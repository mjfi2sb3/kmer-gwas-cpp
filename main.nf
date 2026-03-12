#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KMER_COUNT   } from './modules/kmer_count'
include { MATRIX_MERGE } from './modules/matrix_merge'

// Load utility functions (estimateBins)
evaluate(new File("${projectDir}/lib/utils.nf"))

// ---------------------------------------------------------------------------
// Workflow entry point
// ---------------------------------------------------------------------------
workflow {

    // -- Bin estimation hint (informational only; does not override num_bins) --
    estimateBins(
        params.genome_size,
        params.sequencing_depth,
        params.available_ram_per_node,
        /*bytes_per_kmer=*/ 8,
        params.num_bins
    )

    // -- Stage 1: k-mer counting, one job per accession --
    ch_accessions = Channel
        .fromPath(params.accessions_file)
        .splitText()
        .map { it.trim() }
        .filter { it }                  // skip blank lines

    KMER_COUNT(ch_accessions, params.num_bins)

    // -- Fan-in: wait for ALL accessions, then fire once per bin --
    //
    // Strategy: each finished accession emits num_bins sentinel tuples
    // (bin_idx, kmer_count_root). groupTuple() collects one entry per
    // accession for every bin_idx; the channel only closes (and each bin
    // fires) once all accessions have completed Stage 1.
    //
    def kmer_count_root = "${params.output_dir}/kmer_count"
    def n_bins          = params.num_bins as Integer

    ch_bin_signals = KMER_COUNT.out.accession_dir
        .flatMap { accession, acc_dir ->
            (0..<n_bins).collect { bin_idx ->
                tuple(bin_idx, kmer_count_root)
            }
        }

    ch_bins_ready = ch_bin_signals
        .groupTuple()                           // group by bin_idx; size == num_accessions
        .map { bin_idx, roots -> tuple(bin_idx, roots[0]) }

    // -- Stage 2: matrix merge, one job per bin --
    MATRIX_MERGE(ch_bins_ready, file(params.accessions_file))
}
