process MATRIX_MERGE {
    tag "bin_${bin_idx}"
    publishDir "${params.output_dir}/matrix", mode: 'copy', overwrite: true

    input:
        tuple val(bin_idx), val(kmer_count_root)
        val accessions_file

    output:
        path "matrix_*/${bin_idx}_matrix.tsv.gz", emit: matrix_file
        path "matrix_*/${bin_idx}_core.txt.gz",   emit: core_file, optional: true

    script:
    """
    g++ -std=c++17 -O3 -march=native -pthread -o matrix_merge \
        /opt/kmer-gwas/src/matrix_merge.cpp \
        /opt/kmer-gwas/src/mmap_io.cpp

    mkdir -p extracted/
    for tarball in ${kmer_count_root}/*.tar; do
        acc=\$(basename "\${tarball}" .tar)
        tar -xf "\${tarball}" -C extracted/ "\${acc}/${bin_idx}_nr.bin" 2>/dev/null || true
    done

    ./matrix_merge \\
        --input      extracted/ \\
        --accessions ${accessions_file} \\
        --index      ${bin_idx} \\
        --threshold  ${params.threshold} \\
        --delimiter  ${params.delimiter} \\
        --count      ${params.count} \\
        --core       ${params.core} \\
        --bins       ${params.num_bins} \\
        --threads    ${task.cpus}

    rm -rf extracted/

    if command -v pigz > /dev/null 2>&1; then
        pigz -p ${task.cpus} matrix_*/${bin_idx}_matrix.tsv
        pigz -p ${task.cpus} matrix_*/${bin_idx}_core.txt 2>/dev/null || true
    else
        gzip matrix_*/${bin_idx}_matrix.tsv
        gzip matrix_*/${bin_idx}_core.txt 2>/dev/null || true
    fi
    """
}
