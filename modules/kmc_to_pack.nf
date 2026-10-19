// Stage 1 (kmc input mode) — convert one accession's pre-built KMC database into
// a self-indexed pack, instead of counting k-mers from reads.
//
// The pack is the same `<accession>.kbin` that KMER_COUNT produces (same
// canonical keys, binning and count saturation) and is published to the same
// directory with the same publish mode, so Stage 2 runs unchanged for both input
// modes. Every k-mer in the database is kept unless --min_kmer_count is above
// the database's own -ci (see src/kmc_to_pack.cpp).
//
// The database is read through the KMC C++ API compiled from src/kmc_api/
// (GNU GPL 3), so no KMC installation is needed at run time. k is compiled in
// (-DKMER_K); kmc_to_pack rejects a database built at a different k, so this task
// fails instead of writing a wrong pack.
process KMC_TO_PACK {
    tag "${accession}"

    // Same path as KMER_COUNT: MATRIX_MERGE reads this directory (see main.nf).
    publishDir "${params.output_dir}/kmer_count_k${params.kmer_size}", mode: params.kmer_count_publish_mode, overwrite: true, pattern: "*.kbin"

    input:
        val accession
        val num_bins
        val kmc_dir

    output:
        tuple val(accession), path("${accession}.kbin"), emit: accession_pack
        path "${accession}.pack_stats.tsv",               emit: stats

    script:
    // As in KMER_COUNT, a budget of 0 auto-sizes the accumulation memory from the
    // limit actually enforced on the job; it only affects spilling, never the pack.
    // With the local executor, KMER_GWAS_TASK_MEMORY passes the task's memory
    // request as one more limit, because local tasks share one machine.
    """
    for ext in kmc_pre kmc_suf; do
        [ -f "${kmc_dir}/${accession}.\${ext}" ] || {
            echo "ERROR: missing KMC database file ${kmc_dir}/${accession}.\${ext}" >&2; exit 1; }
    done

    g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=${params.kmer_size} -o kmc_to_pack \\
        ${params.src_dir}/kmc_to_pack.cpp \\
        ${params.src_dir}/kmc_api/kmc_file.cpp \\
        ${params.src_dir}/kmc_api/kmer_api.cpp \\
        ${params.src_dir}/kmc_api/mmer.cpp

    ${task.executor == 'local' && task.memory ? "export KMER_GWAS_TASK_MEMORY=${task.memory.toBytes()}" : ''}
    ./kmc_to_pack ${kmc_dir}/${accession} ${num_bins} ./ ${accession} \\
        ${params.min_kmer_count != null ? params.min_kmer_count : 0} ${params.kmer_count_budget_gb}
    """
}
