// Stage 2 — build one bin's presence/absence matrix by k-way merging the
// matching bin slice from every accession's pack.
//
// Each pack is self-indexed, so matrix_merge seeks directly to the bin it needs
// inside every pack; nothing is unpacked to disk first and no per-bin temporary
// files are created.
//
// The pack directory is passed as a value (an absolute path the task reads
// directly), not a Nextflow `path` input. With large cohorts, staging the packs
// individually would create one symlink per accession for every bin task; and the
// directory is the published output_dir, which is mutated after a run (re-publish
// timestamps, and --compress_kbin_packs rewriting each .kbin to .kbin.gz). Hashing
// it as a `path` broke -resume — see the note in main.nf. The `# cache signature`
// line below folds the Stage-1 pack-affecting param that is not otherwise in this
// script (min_kmer_count) into the task hash, so changing it re-runs the merge.
process MATRIX_MERGE {
    tag "bin_${bin_idx}"
    publishDir "${params.output_dir}/matrix", mode: params.matrix_publish_mode, overwrite: true

    input:
        tuple val(bin_idx), val(kmer_dir)
        val accessions_file

    output:
        path "matrix_*/${bin_idx}_matrix.tsv.gz", emit: matrix_file
        path "matrix_*/${bin_idx}_core.txt.gz",   emit: core_file, optional: true

    script:
    """
    # cache signature (do not remove): folds Stage-1 pack-affecting params into this
    # task's hash so -resume re-runs the merge when they change. min_kmer_count=${params.min_kmer_count}
    #
    # -DKMER_K must match the value KMER_COUNT used. A mismatch is caught at run
    # time by the pack footer, but building consistently avoids the error.
    g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=${params.kmer_size} -o matrix_merge \
        ${params.src_dir}/matrix_merge.cpp

    ./matrix_merge \\
        --input      ${kmer_dir}/ \\
        --accessions ${accessions_file} \\
        --index      ${bin_idx} \\
        --threshold  ${params.threshold} \\
        --delimiter  ${params.delimiter} \\
        --encoding   ${params.encoding} \\
        --count      ${params.keep_kmer_counts} \\
        --core       ${params.write_core_kmers ? 'y' : 'n'} \\
        --bins       ${params.num_bins} \\
        --threads    ${task.cpus}

    if command -v pigz > /dev/null 2>&1; then
        pigz -p ${task.cpus} matrix_*/${bin_idx}_matrix.tsv
        pigz -p ${task.cpus} matrix_*/${bin_idx}_core.txt 2>/dev/null || true
    else
        gzip matrix_*/${bin_idx}_matrix.tsv
        gzip matrix_*/${bin_idx}_core.txt 2>/dev/null || true
    fi
    """
}
