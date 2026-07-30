#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KMER_COUNT   } from './modules/kmer_count'
include { MATRIX_MERGE } from './modules/matrix_merge'
include { BUILD_TOOLS  } from './modules/build_tools'

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    Usage:
        nextflow run main.nf [options]

    Options:
        --accessions_file        Path to file listing accession IDs, one per line   [default: ${params.accessions_file}]
        --data_dir               Directory containing paired FASTQ files            [default: ${params.data_dir}]
        --output_dir             Output directory                                   [default: ${params.output_dir}]
        --kmer_size              k-mer length; must be ODD, 15..63                  [default: ${params.kmer_size}]
                                 k <= 32 uses a compact 10-byte record instead of 18,
                                 cutting Stage 1 output by ~44%. Compiled in, so there
                                 is no runtime cost to changing it.
                                 NOTE: changing k invalidates existing bin files.
        --num_bins               Number of k-mer bins (parallelism / output size)   [default: ${params.num_bins}]
        --min_kmer_count         Drop k-mers seen fewer times in an accession       [default: ${params.min_kmer_count}]
        --threshold              Two-sided MAF filter on accession occurrence      [default: ${params.threshold}]
                                 Keeps k-mers in [threshold, n_acc - threshold] accessions.
                                 Filters on HOW MANY ACCESSIONS carry the k-mer, not on
                                 its count within an accession. 0 disables the filter.
        --count                  'y' = counts, 'n' = presence/absence               [default: ${params.count}]
        --delimiter              Value separator (text encoding): 'tab','space','none' [default: ${params.delimiter}]
                                 'none' concatenates and is presence/absence only.
        --encoding               Matrix encoding: 'text' or 'bits'                  [default: ${params.encoding}]
                                 'bits' packs presence/absence 1 bit per accession
                                 as hex (~8x smaller than tab); requires --count n.
        --core                   'y' = write core k-mers file, 'n' = skip           [default: ${params.core}]
                                 Core k-mers (present in ALL accessions) are excluded
                                 from the matrix — they carry no association signal.
        --matrix_merge_cpus      Number of threads for MATRIX_MERGE                 [default: ${params.matrix_merge_cpus}]
        --kmer_count_memory      RAM requested per KMER_COUNT job                    [default: ${params.kmer_count_memory}]
                                 (a scheduler request; the counting budget auto-sizes
                                 from it — see --kmer_count_budget_gb)
        --kmer_count_budget_gb   Stage 1 accumulation budget in GB; 0 = auto,       [default: ${params.kmer_count_budget_gb}]
                                 a fraction of the RAM actually enforced on the job
        --kmer_count_read_threads  Decompress the two mate files concurrently       [default: ${params.kmer_count_read_threads}]
                                 (set 1 for a single spinning disk)
        --matrix_merge_memory    RAM requested per MATRIX_MERGE job                 [default: ${params.matrix_merge_memory}]
                                 Use dot notation: 64.MB, 120.GB, 256.GB
                                 (NOT '120 GB' — the space form fails on the CLI)
        --kmer_count_time        Wallclock time limit for KMER_COUNT                [default: ${params.kmer_count_time}]
        --matrix_merge_time      Wallclock time limit for MATRIX_MERGE              [default: ${params.matrix_merge_time}]
                                 Use quoted string: '5h', '10h', '1d', '2h 30m'
        --clusterOptions         Extra SLURM options passed to all jobs             [default: none]
                                 Use = syntax: --clusterOptions='--account=myproject --partition=highmem'
                                 To request all node RAM: --clusterOptions='--account=myproject --mem=0'
                                 (--mem=0 takes precedence over --kmer_count_memory/--matrix_merge_memory)
        --singularity_cache_dir  Local path for Singularity image cache             [default: .singularity/]
        --queue_size             Max SLURM jobs submitted (queued+running) at once  [default: ${params.queue_size}]
        --cleanup                Delete the work dir on successful completion        [default: false]
                                 Default false keeps it so -resume can skip finished work;
                                 set --cleanup true to delete on success (disables -resume).
        --publish_mode           How Stage 1 packs reach output_dir:                [default: link]
                                 link (hardlink; enables -resume, same filesystem), copy (~2x storage),
                                 or move (lowest footprint, no -resume).

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
      kmer_size              : ${c_val}${params.kmer_size}${c_reset}
      num_bins               : ${c_val}${params.num_bins}${c_reset}
      min_kmer_count         : ${c_val}${params.min_kmer_count}${c_reset}
      threshold              : ${c_val}${params.threshold}${c_reset}
      count                  : ${c_val}${params.count}${c_reset}
      delimiter              : ${c_val}${params.delimiter}${c_reset}
      encoding               : ${c_val}${params.encoding}${c_reset}
      core                   : ${c_val}${params.core}${c_reset}
      matrix_merge_cpus      : ${c_val}${params.matrix_merge_cpus}${c_reset}
      kmer_count_memory      : ${c_val}${params.kmer_count_memory.toGiga()}.GB${c_reset}
      matrix_merge_memory    : ${c_val}${params.matrix_merge_memory.toGiga()}.GB${c_reset}
      kmer_count_time        : ${c_val}${params.kmer_count_time}${c_reset}
      matrix_merge_time      : ${c_val}${params.matrix_merge_time}${c_reset}
      clusterOptions         : ${c_val}${params.clusterOptions ?: '(none)'}${c_reset}
      singularity_cache_dir  : ${c_val}${params.singularity_cache_dir}${c_reset}
      cleanup                : ${c_val}${params.cleanup}${c_reset}
    ${c_head}Output${c_reset}
      output_dir             : ${c_val}${params.output_dir}${c_reset}
    ${c_banner}-----------------------------------------${c_reset}
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Run manifest
// ---------------------------------------------------------------------------
// Records the parameters a results directory was produced with. k in particular
// is not recoverable from the output files themselves, and reading bin files at
// the wrong k yields silent nonsense rather than an error — so the manifest is
// the record of what a given results/ actually contains.
def writeManifest(String accessions_file, String data_dir, int n_accessions) {
    def dir = file(params.output_dir)
    dir.mkdirs()
    file("${params.output_dir}/run_manifest.txt").text = """\
        # kmer-GWAS run manifest
        # Written at launch; describes how this results directory was produced.
        pipeline_version    = ${workflow.manifest.version ?: 'n/a'}
        run_name            = ${workflow.runName}
        started             = ${workflow.start}
        nextflow_version    = ${workflow.nextflow.version}
        pipeline_revision   = ${workflow.revision ?: 'n/a'}
        commit_id           = ${workflow.commitId ?: 'n/a'}
        command_line        = ${workflow.commandLine}

        kmer_size           = ${params.kmer_size}
        num_bins            = ${params.num_bins}
        min_kmer_count      = ${params.min_kmer_count}
        threshold           = ${params.threshold}
        count               = ${params.count}
        delimiter           = ${params.delimiter}
        encoding            = ${params.encoding}
        core                = ${params.core}

        accessions_file     = ${accessions_file}
        n_accessions        = ${n_accessions}
        data_dir            = ${data_dir}
        output_dir          = ${params.output_dir}

        # Stage 1 output lives in kmer_count_k${params.kmer_size}/ — bin files store
        # fixed-width ${params.kmer_size}-mers and cannot be read at any other k.
        """.stripIndent()
}

// ---------------------------------------------------------------------------
// Completion summary
// ---------------------------------------------------------------------------
workflow.onComplete {
    // NB: work-dir cleanup is handled by the native `cleanup` directive in
    // nextflow.config (cleanup = params.cleanup), which is symlink-safe. Do NOT
    // reintroduce `file(workflow.workDir).deleteDir()` here: MATRIX_MERGE stages
    // the published kmer_count_k<k>/ directory into its work dir as a symlink
    // (path(kmer_dir)), and Groovy's deleteDir() follows that directory symlink
    // — recursing into results/ and deleting the real .kbin packs. That was the
    // cause of the data loss where --cleanup true wiped kmer_count_k*/*.kbin.
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


    // -- Validate k up front -------------------------------------------------
    // k is compiled in as -DKMER_K, so an invalid value would otherwise surface
    // as a static_assert failure inside every SLURM job rather than here.
    // k must be odd (an even k lets a sequence be its own reverse complement,
    // making the canonical form ambiguous) and fit two 64-bit words.
    def kmer_size = params.kmer_size as Integer
    if (kmer_size % 2 == 0 || kmer_size < 15 || kmer_size > 63) {
        exit 1, "ERROR: --kmer_size must be an ODD number between 15 and 63 (got ${params.kmer_size})"
    }

    // Validate the output-format options up front, so a bad combination fails at
    // launch rather than in every MATRIX_MERGE job. Mirrors the C++ checks.
    if (!(params.delimiter in ['tab', 'space', 'none']))
        exit 1, "ERROR: --delimiter must be 'tab', 'space' or 'none' (got '${params.delimiter}'). 'bits' is now --encoding bits."
    if (!(params.encoding in ['text', 'bits']))
        exit 1, "ERROR: --encoding must be 'text' or 'bits' (got '${params.encoding}')."
    if (params.encoding == 'bits' && params.count == 'y')
        exit 1, "ERROR: --encoding bits is presence/absence only; it cannot be combined with --count y."
    if (params.encoding == 'text' && params.delimiter == 'none' && params.count == 'y')
        exit 1, "ERROR: --delimiter none cannot be combined with --count y (multi-digit counts would be unparseable). Use --delimiter tab or space."

    // Validate the publish mode, and for 'link' fail here if output_dir and the
    // work dir are on different filesystems. A hard link cannot cross filesystems,
    // so Nextflow would otherwise fail the FIRST pack publish mid-run with a
    // cryptic error; this stops at launch with the fix (use --publish_mode copy).
    if (!(params.publish_mode in ['link', 'copy', 'move']))
        exit 1, "ERROR: --publish_mode must be 'link', 'copy' or 'move' (got '${params.publish_mode}')."
    if (params.publish_mode == 'link') {
        def outDev = null, workDev = null
        try {
            def devOf = { p ->
                def f = new File(p.toString())
                while (f != null && !f.exists()) f = f.parentFile   // nearest existing ancestor
                f == null ? null : java.nio.file.Files.getAttribute(f.toPath(), 'unix:dev')
            }
            outDev  = devOf(params.output_dir)
            workDev = devOf(workflow.workDir)
        } catch (Exception e) { /* cannot determine filesystems; let Nextflow surface any link error */ }
        if (outDev != null && workDev != null && outDev != workDev)
            exit 1, "ERROR: --publish_mode 'link' needs --output_dir and the work dir on the same filesystem, " +
                    "but they are on different ones (output_dir=${params.output_dir}, work=${workflow.workDir}). " +
                    "Use --publish_mode copy, or point --output_dir at the same filesystem as the work dir."
    }

    // Resolve to absolute paths — relative inputs break inside Singularity
    // containers and SLURM work directories
    def accessions_file = file(params.accessions_file).toAbsolutePath().toString()
    def data_dir        = file(params.data_dir).toAbsolutePath().toString()

    paramSummary(accessions_file, data_dir)

    def n_accessions = file(accessions_file).readLines().findAll { it.trim() }.size()
    writeManifest(accessions_file, data_dir, n_accessions)

    // -- Stage 1: k-mer counting, one job per accession --
    ch_accessions = Channel
        .fromPath(params.accessions_file)
        .splitText()
        .map { it.trim() }
        .filter { it }                  // skip blank lines

    // Build the bits_to_text converter into results/bin/ so the output
    // directory ships with the tool needed to read its bit-packed matrices.
    BUILD_TOOLS()

    KMER_COUNT(ch_accessions, params.num_bins, data_dir)

    // -- Fan-in: wait for ALL accessions, then fire once per bin --
    //
    // Strategy: each finished accession emits num_bins sentinel tuples
    // (bin_idx, kmer_count_root). groupTuple() collects one entry per
    // accession for every bin_idx; the channel only closes (and each bin
    // fires) once all accessions have completed Stage 1.
    //
    // Must match KMER_COUNT's publishDir exactly — k is part of the path so that
    // a results directory from a previous run at a different k cannot be picked
    // up silently.
    def kmer_count_root = "${params.output_dir}/kmer_count_k${params.kmer_size}"
    def n_bins          = params.num_bins as Integer

    ch_bin_signals = KMER_COUNT.out.accession_pack
        .flatMap { accession, pack ->
            (0..<n_bins).collect { bin_idx ->
                tuple(bin_idx, kmer_count_root)
            }
        }

    // One directory path per bin, not one path per accession: MATRIX_MERGE stages
    // it as a single symlink. Passing the packs individually through the channel
    // would cost one symlink per accession per bin task (millions of inodes for a
    // large cohort at queueSize 200), which is worse than the extraction it replaced.
    ch_bins_ready = ch_bin_signals
        .groupTuple()                           // group by bin_idx; size == num_accessions
        .map { bin_idx, roots -> tuple(bin_idx, file(roots[0])) }

    // -- Stage 2: matrix merge, one job per bin --
    MATRIX_MERGE(ch_bins_ready, file(accessions_file))
}
