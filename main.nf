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
def paramSummary(String accessions_file, String data_dir) {
    def c_reset  = "\033[0m"
    def c_banner = "\033[1;36m"   // bold cyan  — banner & dividers
    def c_head   = "\033[1;33m"   // bold yellow — section headers
    def c_val    = "\033[0;32m"   // green       — values

    log.info """
    ${c_banner}=========================================${c_reset}
    ${c_banner} k m e r - G W A S  p i p e l i n e${c_reset}
    ${c_banner}=========================================${c_reset}
    ${c_head}Input${c_reset}
      accessions_file        : ${c_val}${accessions_file}${c_reset}
      data_dir               : ${c_val}${data_dir}${c_reset}
    ${c_head}Pipeline${c_reset}
      num_bins               : ${c_val}${params.num_bins}${c_reset}
      threshold              : ${c_val}${params.threshold}${c_reset}
      count                  : ${c_val}${params.count}${c_reset}
      delimiter              : ${c_val}${params.delimiter}${c_reset}
      core                   : ${c_val}${params.core}${c_reset}
      matrix_merge_cpus      : ${c_val}${params.matrix_merge_cpus}${c_reset}
    ${c_head}Output${c_reset}
      output_dir             : ${c_val}${params.output_dir}${c_reset}
    ${c_banner}-----------------------------------------${c_reset}
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Completion summary
// ---------------------------------------------------------------------------
workflow.onComplete {
    def c_reset  = "\033[0m"
    def c_banner = workflow.success ? "\033[1;36m" : "\033[1;31m"  // bold cyan or bold red
    def c_head   = "\033[1;33m"                                     // bold yellow
    def c_val    = "\033[0;32m"                                     // green
    def c_fail   = "\033[0;31m"                                     // red

    log.info """
    ${c_banner}=========================================${c_reset}
    ${c_banner} Pipeline ${workflow.success ? 'completed' : 'FAILED'}${c_reset}
    ${c_banner}=========================================${c_reset}
    ${c_head}Run${c_reset}
      Run name   : ${c_val}${workflow.runName}${c_reset}
      Completed  : ${c_val}${workflow.complete}${c_reset}
      Duration   : ${c_val}${workflow.duration}${c_reset}
      CPU hours  : ${c_val}${workflow.stats.computeTimeFmt}${c_reset}
    ${c_head}Tasks${c_reset}
      Succeeded  : ${c_val}${workflow.stats.succeededCount}${c_reset}
      Cached     : ${c_val}${workflow.stats.cachedCount}${c_reset}
      Failed     : ${workflow.stats.failedCount > 0 ? c_fail : c_val}${workflow.stats.failedCount}${c_reset}
    ${c_head}Output${c_reset}
      output_dir : ${c_val}${params.output_dir}${c_reset}
    ${c_banner}-----------------------------------------${c_reset}
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

    // Resolve to absolute paths — relative inputs break inside Singularity
    // containers and SLURM work directories
    def accessions_file = file(params.accessions_file).toAbsolutePath().toString()
    def data_dir        = file(params.data_dir).toAbsolutePath().toString()

    paramSummary(accessions_file, data_dir)

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
    MATRIX_MERGE(ch_bins_ready, file(accessions_file))
}
