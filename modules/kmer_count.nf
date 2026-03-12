process KMER_COUNT {
    tag "${accession}"
    publishDir "${params.output_dir}/kmer_count", mode: 'copy', overwrite: true

    input:
        val accession
        val num_bins

    output:
        tuple val(accession), path("${accession}"), emit: accession_dir

    script:
    """
    mkdir -p data
    ln -sf ${params.data_dir}/${accession}_1.fq data/${accession}_1.fq
    ln -sf ${params.data_dir}/${accession}_2.fq data/${accession}_2.fq
    ${projectDir}/build/kmer_count ${accession} ${num_bins} ./
    """
}
