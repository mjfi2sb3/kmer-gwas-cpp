// Optional archival step: gzip the published Stage 1 packs, one job per accession.
//
// Runs ONLY after Stage 2 (see main.nf): the process is gated on every
// MATRIX_MERGE task, and Stage 2 already depends on all of Stage 1, so no pack is
// ever compressed while it might still be read. If any Stage 1/2 task fails, the
// run aborts before this ever launches, leaving the packs uncompressed.
//
// The gain is modest (~1.5-2x, varying with k and coverage): the 2-bit packed
// keys are near-incompressible, so most of it comes from key padding bits and the
// small count column. It rewrites <accession>.kbin -> <accession>.kbin.gz in the PUBLISHED pack
// directory (output_dir/kmer_count_k<k>/), in place with pigz. Note this cannot
// preserve a hard link: pigz writes a new .gz and removes the .kbin, so if the
// uncompressed pack is also hard-linked in the work dir (kmer_count_publish_mode 'link' with
// cleanup=false) the data is not freed until the work dir is removed. main.nf
// warns about that at launch. The read path (Stage 2, kbin_dump) always expects
// uncompressed packs, so re-running Stage 2 means decompressing first.
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
