/**
 * Estimate the recommended number of k-mer bins based on genome and sequencing parameters.
 *
 * Formula (rough):
 *   total_kmers = genome_size * sequencing_depth
 *   ram_needed_gb = total_kmers * num_accessions * bytes_per_kmer / 1024^3
 *   recommended_bins = ceil(ram_needed_gb / available_ram_gb)
 *
 * Note: does NOT account for k-mer deduplication (~60-70% for plant genomes).
 * Future work: add --auto_bins flag to override num_bins automatically.
 */
def estimateBins(genome_size, sequencing_depth, available_ram_gb, bytes_per_kmer, user_bins) {
    if (!genome_size || !sequencing_depth || !available_ram_gb) {
        log.info "[estimateBins] Skipping — genome_size, sequencing_depth, or available_ram_per_node not set. Using num_bins=${user_bins}"
        return
    }
    def num_accessions = file(params.accessions_file).readLines().findAll { it.trim() }.size()
    long   total_kmers  = (genome_size as Long) * (sequencing_depth as Integer)
    double ram_gb       = total_kmers * num_accessions * (bytes_per_kmer as Integer) / (1024.0**3)
    long   recommended  = Math.ceil(ram_gb / (available_ram_gb as Double)) as Long

    log.info """
    [estimateBins] genome=${genome_size}bp depth=${sequencing_depth}x accessions=${num_accessions}
      => estimated RAM without binning: ${String.format('%.1f', ram_gb)} GB
      => recommended num_bins: ${recommended}  (you specified: ${user_bins})
    """.stripIndent()
}
