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
    def extract = params.large_cohort ? """
        mkdir -p extracted/
        for tarball in ${kmer_count_root}/*.tar; do
            acc=\$(basename "\${tarball}" .tar)
            tar -xf "\${tarball}" "\${acc}/${bin_idx}_nr.bin" -C extracted/ 2>/dev/null || true
        done
    """ : ""
    def input_dir = params.large_cohort ? "extracted/" : "${kmer_count_root}/"
    def cleanup   = params.large_cohort ? "rm -rf extracted/" : ""
    """
    g++ -std=c++17 -O3 -march=native -pthread -o matrix_merge \
        /opt/kmer-gwas/src/matrix_merge.cpp \
        /opt/kmer-gwas/src/mmap_io.cpp

    ${extract}

    ./matrix_merge \\
        --input      ${input_dir} \\
        --accessions ${accessions_file} \\
        --index      ${bin_idx} \\
        --threshold  ${params.threshold} \\
        --delimiter  ${params.delimiter} \\
        --count      ${params.count} \\
        --core       ${params.core} \\
        --bins       ${params.num_bins} \\
        --threads    ${task.cpus}

    ${cleanup}
    """
}
