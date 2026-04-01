#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KMER_COUNT   } from './modules/kmer_count'
include { MATRIX_MERGE } from './modules/matrix_merge'

// Load utility functions (estimateBins)
//evaluate(new File("${projectDir}/lib/utils.nf"))

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    Usage:
        nextflow run main.nf [options]

    Options:
        --accessions_file   Path to file listing accession IDs, one per line  [default: ${params.accessions_file}]
        --data_dir          Directory containing paired FASTQ files            [default: ${params.data_dir}]
        --output_dir        Output directory                                   [default: ${params.output_dir}]
        --num_bins          Number of k-mer bins                               [default: ${params.num_bins}]
        --threshold         Minimum count threshold for k-mer inclusion        [default: ${params.threshold}]
        --count             'y' = counts, 'n' = presence/absence               [default: ${params.count}]
        --delimiter         Column delimiter: 'tab' or 'none'                  [default: ${params.delimiter}]
        --core              'y' = write core k-mers file, 'n' = skip           [default: ${params.core}]
        --matrix_merge_cpus Number of threads for matrix_merge                 [default: ${params.matrix_merge_cpus}]

    Profiles:
        -profile standard           Run locally
        -profile slurm              Submit jobs to SLURM (module-system environment)
        -profile slurm_container    Submit jobs to SLURM using Singularity container (recommended)

    Example:
        nextflow run main.nf -profile slurm_container --accessions_file samples.txt --data_dir /path/to/fastq
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


// ---------------------------------------------------------------------------
// Parameter summary
// ---------------------------------------------------------------------------
def paramSummary(String data_dir) {
    log.info """
    =========================================
     k m e r - G W A S  p i p e l i n e
    =========================================
    Input
      accessions_file        : ${params.accessions_file}
      data_dir               : ${data_dir}
    Pipeline
      num_bins               : ${params.num_bins}
      threshold              : ${params.threshold}
      count                  : ${params.count}
      delimiter              : ${params.delimiter}
      core                   : ${params.core}
      matrix_merge_cpus      : ${params.matrix_merge_cpus}
    Output
      output_dir             : ${params.output_dir}
    -----------------------------------------
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Workflow entry point
// ---------------------------------------------------------------------------
workflow {

    // -- Bin estimation hint (informational only; does not override num_bins) --
    //estimateBins(
    //    params.genome_size,
    //    params.sequencing_depth,
    //    params.available_ram_per_node,
    //    /*bytes_per_kmer=*/ 8,
    //    params.num_bins
    //)

    // Resolve data_dir to absolute path — relative inputs (e.g. ./data) break
    // inside Singularity containers and SLURM work directories
    def data_dir = file(params.data_dir).toAbsolutePath().toString()

    paramSummary(data_dir)

    // -- Stage 1: k-mer counting, one job per accession --
    ch_accessions = Channel
        .fromPath(params.accessions_file)
        .splitText()
        .map { it.trim() }
        .filter { it }                  // skip blank lines

    KMER_COUNT(ch_accessions, params.num_bins, data_dir)

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
