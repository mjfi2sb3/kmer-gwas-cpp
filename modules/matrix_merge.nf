process MATRIX_MERGE {
    tag "bin_${bin_idx}"
    publishDir "${params.output_dir}/matrix", mode: 'move', overwrite: true

    input:
        tuple val(bin_idx), val(kmer_count_root)
        val accessions_file

    output:
        path "matrix_*/${bin_idx}_matrix.tsv", emit: matrix_file
        path "matrix_*/${bin_idx}_core.txt",   emit: core_file, optional: true

    script:
    """
    # Load gcc if not available or below minimum version (9)
    if ! command -v g++ &>/dev/null || [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        module load gcc/12.2.0 2>/dev/null || true
    fi
    if ! command -v g++ &>/dev/null; then
        echo "WARNING: g++ not available after module load attempt" >&2
    elif [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        echo "WARNING: g++ \$(g++ -dumpversion) < 9 -- compilation may fail" >&2
    fi

    g++ -std=c++17 -O3 -march=native -pthread -o matrix_merge \
        ${projectDir}/src/matrix_merge.cpp \
        ${projectDir}/src/mmap_io.cpp

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
