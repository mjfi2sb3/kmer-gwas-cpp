// Stage 2 — build one bin's matrix by k-way merging every accession's pack slice.
//
// The previous version extracted `<acc>/<bin>_nr.bin` out of each accession's
// tar into a local `extracted/` directory: 12,600 directories plus 12,600 files
// PER bin task, ~5M live inodes at queueSize 200. That step is gone entirely —
// matrix_merge seeks straight to the bin it needs inside each pack.
//
// `kmer_dir` is staged as ONE directory symlink rather than 12,600 individual
// file inputs. Staging the files individually would cost 12,600 symlinks per
// bin task, which is worse than what it replaced.
process MATRIX_MERGE {
    tag "bin_${bin_idx}"
    publishDir "${params.output_dir}/matrix", mode: 'copy', overwrite: true

    input:
        tuple val(bin_idx), path(kmer_dir)
        val accessions_file

    output:
        path "matrix_*/${bin_idx}_matrix.tsv.gz", emit: matrix_file
        path "matrix_*/${bin_idx}_core.txt.gz",   emit: core_file, optional: true

    script:
    """
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
        --count      ${params.count} \\
        --core       ${params.core} \\
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
