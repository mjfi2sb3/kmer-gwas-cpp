// Optional archival step: gzip the published Stage 1 packs after both stages
// complete. Compression is disabled if any Stage 1/2 task fails.
//
// With --kmer_count_publish_mode link and --cleanup false, the uncompressed packs
// remain hard-linked in the work dir, so space is only reclaimed after the work
// dir is removed. See main.nf for the launch warning.
//
// Re-running Stage 2 requires decompressing the packs first.
process COMPRESS_PACKS {
    tag "${accession}"

    input:
        tuple val(accession), val(gate)   // gate: all MATRIX_MERGE outputs (collected)
        val   kmer_count_dir

    script:
    """
    pack="${kmer_count_dir}/${accession}.kbin"
    if [ -f "\$pack" ]; then
        echo "COMPRESS_PACKS: pigz -p ${task.cpus} \$pack"
        pigz -p ${task.cpus} -f "\$pack"
    else
        echo "COMPRESS_PACKS: \$pack not found (already compressed or moved) - skipping"
    fi
    """
}
