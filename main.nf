#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KMER_COUNT    } from './modules/kmer_count'
include { MATRIX_MERGE  } from './modules/matrix_merge'
include { BUILD_TOOLS   } from './modules/build_tools'
include { COMPRESS_PACKS } from './modules/compress_packs'

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    Usage:
        nextflow run main.nf [options]

    Options:

      Input / output
        --accessions_file        Path to file listing accession IDs, one per line   [default: ${params.accessions_file}]
        --data_dir               Directory containing paired FASTQ files            [default: ${params.data_dir}]
        --output_dir             Output directory                                   [default: ${params.output_dir}]

      k-mers and matrix content
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
        --keep_kmer_counts       'y' = keep actual kmer counts, 'n' = presence/absence               [default: ${params.count}]
        --delimiter              Value separator (text encoding): 'tab','space','none' [default: ${params.delimiter}]
                                 'none' concatenates and is presence/absence only.
        --encoding               Matrix encoding: 'text' or 'bits'                  [default: ${params.encoding}]
                                 'bits' packs presence/absence 1 bit per accession
                                 as hex (~8x smaller than tab); requires --count n.
        --write_core_kmers       Write k-mers present in ALL accessions to a         [default: ${params.write_core_kmers}]
                                 per-bin <bin>_core.txt file. Independent of every
                                 other flag: it only adds this file and never changes
                                 the matrix. false (default) writes no core file.

      Stage 1 resources (KMER_COUNT, one job per accession)
        --kmer_count_memory      RAM requested per KMER_COUNT job                    [default: ${params.kmer_count_memory}]
                                 A scheduler request; the counting budget auto-sizes
                                 from it (see --kmer_count_budget_gb).
        --kmer_count_budget_gb   Stage 1 accumulation budget in GB; 0 = auto,       [default: ${params.kmer_count_budget_gb}]
                                 a fraction of the RAM actually enforced on the job.
        --kmer_count_read_threads  Decompress the two mate files concurrently       [default: ${params.kmer_count_read_threads}]
                                 (set 1 for a single spinning disk).
        --kmer_count_time        Wallclock time limit for KMER_COUNT                [default: ${params.kmer_count_time}]
                                 Use a quoted string: '5h', '10h', '1d', '2h 30m'.

      Stage 2 resources (MATRIX_MERGE, one job per bin)
        --matrix_merge_cpus      CPUs requested per MATRIX_MERGE job                 [default: ${params.matrix_merge_cpus}]
                                 The merge runs in parallel (byte-identical output)
                                 and these cores also feed pigz; returns diminish
                                 past ~8-12. Higher mainly costs RAM.
        --matrix_merge_memory    RAM requested per MATRIX_MERGE job                 [default: ${params.matrix_merge_memory}]
                                 Use dot notation: 64.MB, 120.GB, 256.GB
                                 (NOT '120 GB'; the space form fails on the CLI).
        --matrix_merge_time      Wallclock time limit for MATRIX_MERGE              [default: ${params.matrix_merge_time}]

      Scheduler and execution
        --clusterOptions         Extra SLURM options passed to all jobs             [default: none]
                                 Use = syntax: --clusterOptions='--account=myproject --partition=highmem'
                                 To request all node RAM: --clusterOptions='--account=myproject --mem=0'
                                 (--mem=0 takes precedence over --kmer_count_memory/--matrix_merge_memory)
        --queue_size             Max SLURM jobs submitted (queued+running) at once  [default: ${params.queue_size}]
        --max_retries            Retries for an OOM/timeout-killed job (0 = none)   [default: ${params.max_retries}]
                                 Each retry gets attempt x the base memory and time.
        --max_memory             Memory ceiling any request is clamped to           [default: ${params.max_memory}]
        --max_cpus               CPU ceiling any request is clamped to              [default: ${params.max_cpus}]
        --max_time               Wallclock ceiling any request is clamped to        [default: ${params.max_time}]
                                 Set --max_memory/--max_cpus/--max_time to your largest node.
        --singularity_cache_dir  Local path for Singularity image cache             [default: .singularity/]

      Work directory and publishing
        --cleanup                Delete the work dir on successful completion        [default: false]
                                 Default false keeps it so -resume can skip finished work;
                                 set --cleanup true to delete on success (disables -resume).
        --kmer_count_publish_mode  How Stage 1 packs reach output_dir:               [default: link]
                                 link (hardlink; enables -resume, same filesystem), copy (~2x storage),
                                 or move (lowest footprint, no -resume). (--publish_mode is a
                                 deprecated alias.)
        --matrix_publish_mode    How the Stage 2 matrix reaches output_dir:          [default: link]
                                 same choices as --kmer_count_publish_mode. link is safe here (the
                                 matrix is already gzipped, never re-compressed) and
                                 avoids copying a large result; use copy for an
                                 independent archival copy.
        --compress_kbin_packs    gzip the published .kbin packs after BOTH stages     [default: ${params.compress_kbin_packs}]
                                 finish (modest, ~1.5-2x: the packed keys barely
                                 compress; the gain varies with k and coverage). Runs
                                 only on full success; re-running Stage 2 then needs
                                 them decompressed first. With kmer_count_publish_mode link +
                                 cleanup false it saves nothing until work is removed.

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

// Backwards compatibility: --core was renamed to --write_core_kmers in v3.7.2.
// The value is mapped through in nextflow.config (params are immutable here);
// this just warns so existing commands keep working. Drop the shim later.
if (params.core != null) {
    log.warn "--core is deprecated and was renamed to --write_core_kmers; " +
             "honouring it as --write_core_kmers ${params.write_core_kmers}. " +
             "Please switch, as --core will be removed in a future release."
}
if (params.publish_mode != null) {
    log.warn "--publish_mode is deprecated and was renamed to --kmer_count_publish_mode; " +
             "honouring it as --kmer_count_publish_mode ${params.kmer_count_publish_mode}. " +
             "Please switch, as --publish_mode will be removed in a future release."
}


// ---------------------------------------------------------------------------
// Parameter summary
// ---------------------------------------------------------------------------
def paramSummary(String accessions_file, String data_dir) {
    def c_reset  = "\033[0m"
    def c_banner = "\033[1;36m"   // bold cyan   — banner & dividers
    def c_head   = "\033[1;33m"   // bold yellow — section headers
    def c_val    = "\033[0;32m"   // green       — values

    // Fixed-width label column (fits the longest key, kmer_count_read_threads),
    // so the colons line up without hand-counting spaces per line.
    def W    = 24
    def row  = { k, v -> "      ${k.padRight(W)}: ${c_val}${v}${c_reset}" }
    def head = { t     -> "    ${c_head}${t}${c_reset}" }
    def rule = "    ${c_banner}=========================================${c_reset}"
    def budget = params.kmer_count_budget_gb == 0 ? '0 (auto)' : "${params.kmer_count_budget_gb} GB"

    // Grouped to mirror the --help sections. workflow.* are runtime facts (the
    // actual version/profile/work dir) that params alone cannot show.
    def lines = [
        "",
        rule,
        "    ${c_banner} k m e r - G W A S  p i p e l i n e${c_reset}",
        rule,
        head('Run'),
        row('version',  workflow.manifest.version),
        row('profile',  workflow.profile),
        row('work_dir', workflow.workDir),
        head('Input / output'),
        row('accessions_file', accessions_file),
        row('data_dir',        data_dir),
        row('output_dir',      params.output_dir),
        head('k-mers and matrix content'),
        row('kmer_size',      params.kmer_size),
        row('num_bins',       params.num_bins),
        row('min_kmer_count', params.min_kmer_count),
        row('threshold',      params.threshold),
        row('keep kmer counts',          params.count),
        row('delimiter',      params.delimiter),
        row('encoding',       params.encoding),
        row('write_core_kmers', params.write_core_kmers),
        head('Stage 1 resources (KMER_COUNT)'),
        row('kmer_count_memory',       "${params.kmer_count_memory.toGiga()}.GB"),
        row('kmer_count_budget_gb',    budget),
        row('kmer_count_read_threads', params.kmer_count_read_threads),
        row('kmer_count_time',         params.kmer_count_time),
        head('Stage 2 resources (MATRIX_MERGE)'),
        row('matrix_merge_cpus',   params.matrix_merge_cpus),
        row('matrix_merge_memory', "${params.matrix_merge_memory.toGiga()}.GB"),
        row('matrix_merge_time',   params.matrix_merge_time),
        head('Scheduler and execution'),
        row('clusterOptions',        params.clusterOptions ?: '(none)'),
        row('queue_size',            params.queue_size),
        row('max_retries',           params.max_retries),
        row('max_memory/cpus/time',  "${params.max_memory} / ${params.max_cpus} / ${params.max_time}"),
        row('singularity_cache_dir', params.singularity_cache_dir),
        head('Work directory and publishing'),
        row('kmer_count_publish_mode', params.kmer_count_publish_mode),
        row('matrix_publish_mode', params.matrix_publish_mode),
        row('cleanup',           params.cleanup),
        row('compress_kbin_packs', params.compress_kbin_packs),
        "    ${c_banner}-----------------------------------------${c_reset}",
        ""
    ]
    log.info lines.join('\n')
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
        keep_kmer_counts    = ${params.count}
        delimiter           = ${params.delimiter}
        encoding            = ${params.encoding}
        write_core_kmers    = ${params.write_core_kmers}

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
    // cryptic error; this stops at launch with the fix (use copy).
    if (!(params.kmer_count_publish_mode in ['link', 'copy', 'move']))
        exit 1, "ERROR: --kmer_count_publish_mode must be 'link', 'copy' or 'move' (got '${params.kmer_count_publish_mode}')."
    if (!(params.matrix_publish_mode in ['link', 'copy', 'move']))
        exit 1, "ERROR: --matrix_publish_mode must be 'link', 'copy' or 'move' (got '${params.matrix_publish_mode}')."
    if (params.kmer_count_publish_mode == 'link' || params.matrix_publish_mode == 'link') {
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
            exit 1, "ERROR: 'link' publish mode needs --output_dir and the work dir on the same filesystem, " +
                    "but they are on different ones (output_dir=${params.output_dir}, work=${workflow.workDir}). " +
                    "Set the offending mode (--kmer_count_publish_mode / --matrix_publish_mode) to copy, " +
                    "or point --output_dir at the same filesystem as the work dir."
    }

    // --compress_kbin_packs only reclaims space if the uncompressed packs are
    // actually gone. With kmer_count_publish_mode 'link' + cleanup=false the pack survives as
    // a hard link in the work dir, so compressing the published copy only ADDS the
    // .gz until the work dir is removed. Warn rather than fail (option (a)).
    if (params.compress_kbin_packs && params.kmer_count_publish_mode == 'link' && !params.cleanup)
        log.warn "--compress_kbin_packs is on, but with --kmer_count_publish_mode link and --cleanup false " +
                 "the uncompressed packs remain hard-linked in the work dir, so no space is " +
                 "reclaimed until it is removed. Use --cleanup true (or --kmer_count_publish_mode move)."

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

    // -- Optional: compress the published packs, one job per accession, only AFTER
    //    every bin's MATRIX_MERGE has finished (gated on the collected outputs), so
    //    a pack is never compressed while it might still be read. See the module.
    if (params.compress_kbin_packs) {
        ch_gate = MATRIX_MERGE.out.matrix_file.collect().map { true }
        ch_accs = KMER_COUNT.out.accession_pack.map { acc, pack -> acc }
        COMPRESS_PACKS(ch_accs.combine(ch_gate), kmer_count_root)
    }
}
