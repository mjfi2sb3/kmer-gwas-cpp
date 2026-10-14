// kmc input mode — read the header of every accession's KMC database before any
// conversion starts, write them to kmc_databases.tsv, and stop the run if the
// databases cannot be combined in one matrix: a missing or unreadable database,
// a k-mer length other than --kmer_size, or different -ci, -cx or (when counts
// are kept) counter sizes. --kmc_allow_mixed turns the last three into
// warnings. A --min_kmer_count at least as large as every database's -ci makes
// different -ci values consistent, so they are then allowed.
//
// One job reads all headers; each read touches only the ends of the files.
process KMC_HEADERS {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true, pattern: "kmc_databases.tsv"

    input:
        path accessions_file
        val  kmc_dir

    output:
        path "kmc_databases.tsv", emit: table
        path "kmc_check.txt",     emit: notes

    script:
    """
    g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=${params.kmer_size} -o kmc_to_pack \\
        ${params.src_dir}/kmc_to_pack.cpp \\
        ${params.src_dir}/kmc_api/kmc_file.cpp \\
        ${params.src_dir}/kmc_api/kmer_api.cpp \\
        ${params.src_dir}/kmc_api/mmer.cpp

    ./kmc_to_pack --headers ${accessions_file} ${kmc_dir} kmc_databases.tsv \\
        ${params.min_kmer_count != null ? params.min_kmer_count : 0} \\
        ${params.keep_kmer_counts} ${params.kmc_allow_mixed ? 'y' : 'n'} > kmc_check.txt
    """
}
