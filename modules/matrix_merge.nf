process MATRIX_MERGE {
    tag "bin_${bin_idx}"
    publishDir "${params.output_dir}/matrix", mode: 'copy', overwrite: true

    input:
        tuple val(bin_idx), val(kmer_count_root)
        val accessions_file

    output:
        path "matrix_*/${bin_idx}_matrix.tsv", emit: matrix_file
        path "matrix_*/${bin_idx}_core.txt",   emit: core_file, optional: true

    script:
    """
    g++ -std=c++17 -O3 -march=native -pthread -o matrix_merge \
        /opt/kmer-gwas/src/matrix_merge.cpp \
        /opt/kmer-gwas/src/mmap_io.cpp

    ./matrix_merge \\
        --input      ${kmer_count_root}/ \\
        --accessions ${accessions_file} \\
        --index      ${bin_idx} \\
        --threshold  ${params.threshold} \\
        --delimiter  ${params.delimiter} \\
        --count      ${params.count} \\
        --core       ${params.core} \\
        --bins       ${params.num_bins} \\
        --threads    ${task.cpus}
    """
}
