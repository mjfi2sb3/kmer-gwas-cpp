#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KMER_COUNT    } from './modules/kmer_count'
include { KMC_TO_PACK   } from './modules/kmc_to_pack'
include { KMC_HEADERS   } from './modules/kmc_headers'
include { MATRIX_MERGE  } from './modules/matrix_merge'
include { BUILD_TOOLS   } from './modules/build_tools'
include { COMPRESS_PACKS } from './modules/compress_packs'
include { validateParameters } from 'plugin/nf-schema'

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
        --input_mode             'reads' (paired FASTQ/FASTA) or 'kmc' (pre-built   [default: ${params.input_mode}]
                                 KMC databases, <kmc_dir>/<accession>.kmc_pre/.kmc_suf)
        --kmc_dir                Directory holding the KMC databases (kmc mode)     [default: data_dir]
        --kmc_allow_mixed        kmc mode: warn instead of stopping when databases  [default: ${params.kmc_allow_mixed}]
                                 were built with different -ci, -cx or counter sizes

      k-mers and matrix content
        --kmer_size              k-mer length; must be ODD, 15..63                  [required in reads mode]
                                 In kmc mode it is read from the databases; if you
                                 also pass it, it must match them.
                                 k <= 32 uses a compact 10-byte record instead of 18,
                                 cutting Stage 1 output by ~44%. Compiled in, so there
                                 is no runtime cost to changing it.
                                 NOTE: changing k invalidates existing bin files.
                                 Output k-mers are canonical: of a k-mer and its reverse
                                 complement, the smaller one (A<C<G<T order) is kept.
        --num_bins               Number of k-mer bins (parallelism / output size)   [default: ${params.num_bins}]
        --min_kmer_count         Drop k-mers seen fewer times in an accession       [default: 2 in reads mode]
                                 In kmc mode the databases are kept as they are (their own
                                 -ci applies); a value above their -ci is enforced as well.
        --threshold              Two-sided MAF filter on accession occurrence      [default: ${params.threshold}]
                                 Keeps k-mers in [threshold, n_acc - threshold] accessions.
                                 Filters on HOW MANY ACCESSIONS carry the k-mer, not on
                                 its count within an accession. 0 disables the filter.
        --keep_kmer_counts       'y' = keep raw kmer counts, 'n' = presence/absence   [default: ${params.keep_kmer_counts}]
                                 (--count is a deprecated alias)
        --skip_accession_check   Allow packs whose recorded accession differs from  [default: ${params.skip_accession_check}]
                                 their file name (packs renamed on purpose)
        --delimiter              Value separator (text encoding): 'tab','space','none' [default: ${params.delimiter}]
                                 'none' concatenates and is presence/absence only.
        --encoding               Matrix encoding: 'text' or 'bits'                  [default: ${params.encoding}]
                                 'bits' packs presence/absence 1 bit per accession
                                 as hex (~8x smaller than tab); requires --keep_kmer_counts n.
        --write_core_kmers       Write k-mers present in ALL accessions to a         [default: ${params.write_core_kmers}]
                                 per-bin <bin>_core.txt file. Independent of every
                                 other flag: it only adds this file and never changes
                                 the matrix. false (default) writes no core file.

      Stage 1 resources (KMER_COUNT, one job per accession)
        --kmer_count_cpus        CPUs requested per KMER_COUNT job                   [default: ${params.kmer_count_cpus}]
                                 kmer_count runs this many worker threads.
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
        --singularity_cache_dir  Singularity/Apptainer image cache directory        [default: .singularity/]

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

    Profiles: one executor, optionally with one container engine (comma-separated)
        executor:  local          Run on this machine (the default)
                   slurm          Submit each task as a SLURM job
        container: singularity, apptainer, docker, podman
                   Without one, the host must provide GCC 12+, zlib and pigz.
        Examples:  -profile slurm,singularity   -profile local,docker

    Example:
        nextflow run main.nf -profile slurm,singularity --accessions_file samples.txt --kmer_size 51 --data_dir /path/to/fastq
    """.stripIndent()
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

    // Grouped to mirror the --help sections. workflow.* are runtime facts (the
    // actual version/profile/work dir) that params alone cannot show.
    def lines = [
        "",
        rule,
        "    ${c_banner} k m e r - G W A S  p i p e l i n e${c_reset}",
        rule,
        head.call('Run'),
        row.call('version',  workflow.manifest.version),
        row.call('profile',  workflow.profile),
        row.call('executor', executorName()),
        row.call('container_engine', engineName()),
        row.call('container', containerText()),
        row.call('work_dir', workflow.workDir),
        head.call('Input / output'),
        row.call('accessions_file', accessions_file),
        row.call('data_dir',        params.input_mode == 'kmc' ? '(n/a: kmc mode)' : data_dir),
        row.call('output_dir',      params.output_dir),
        row.call('input_mode',      params.input_mode),
        row.call('kmc_dir',         params.input_mode == 'kmc' ? file(params.kmc_dir ?: data_dir).toAbsolutePath().toString() : '(n/a: reads mode)'),
        head.call('k-mers and matrix content'),
        row.call('kmer_size',      kmerSizeText()),
        row.call('num_bins',       params.num_bins),
        row.call('min_kmer_count', minCountText()),
        row.call('threshold',      params.threshold),
        row.call('keep kmer counts',          params.keep_kmer_counts),
        row.call('delimiter',      params.delimiter),
        row.call('encoding',       params.encoding),
        row.call('write_core_kmers', params.write_core_kmers),
        head.call(params.input_mode == 'kmc' ? 'Stage 1 resources (KMC_TO_PACK)' : 'Stage 1 resources (KMER_COUNT)')
    ] + stage1Rows().collect { k, v -> row.call(k, v) } + [
        head.call('Stage 2 resources (MATRIX_MERGE, COMPRESS_PACKS)')
    ] + stage2Rows().collect { k, v -> row.call(k, v) } + [
        head.call('Scheduler and execution')
    ] + executionRows().collect { k, v -> row.call(k, v) } + [
        head.call('Work directory and publishing'),
        row.call('kmer_count_publish_mode', params.kmer_count_publish_mode),
        row.call('matrix_publish_mode', params.matrix_publish_mode),
        row.call('cleanup',           params.cleanup),
        row.call('compress_kbin_packs', params.compress_kbin_packs),
        "    ${c_banner}-----------------------------------------${c_reset}",
        ""
    ]
    log.info lines.join('\n')
}

// The Stage 1 count filter as it applies in this run.
def minCountText() {
    if (params.input_mode == 'kmc')
        return params.min_kmer_count == null ? "not set: each database's own -ci applies (see kmc_databases.tsv)"
                                             : "${params.min_kmer_count} (applied where it is above a database's -ci)"
    return params.min_kmer_count == null ? '2' : "${params.min_kmer_count}"
}

// The k-mer size, and in kmc mode where it came from.
def kmerSizeText() {
    if (params.input_mode != 'kmc') return "${params.kmer_size}"
    if (params.kmc_kmer_size_from_db) return "${params.kmer_size} (read from the KMC databases)"
    if (params.kmc_detected_kmer_size != null) return "${params.kmer_size} (given; matches the KMC databases)"
    return "${params.kmer_size} (given; not read from the KMC databases)"
}

// One MD5 over the code a run uses: main.nf, nextflow.config, modules/*.nf and the
// C++ sources in src_dir (.cpp, .hpp, .h, Makefile; built binaries are left out).
// Files are taken in sorted order and each file's relative path is part of the
// digest, so two copies of the code give the same value only if these files match.
def codeChecksum() {
    try {
        def md    = java.security.MessageDigest.getInstance('MD5')
        def files = [:]
        def add   = { String tag, File root, Closure keep ->
            if (!root.isDirectory()) return
            root.eachFileRecurse(groovy.io.FileType.FILES) { f ->
                if (keep.call(f)) files["${tag}/${root.toPath().relativize(f.toPath())}".toString()] = f
            }
        }
        def proj = projectDir.toFile()
        ['main.nf', 'nextflow.config', 'nextflow_schema.json'].each { n -> def f = new File(proj, n); if (f.isFile()) files[n] = f }
        add.call('modules', new File(proj, 'modules'), { File f -> f.name.endsWith('.nf') })
        add.call('conf',    new File(proj, 'conf'),    { File f -> f.name.endsWith('.config') })
        add.call('src', new File(params.src_dir.toString()),
            { File f -> f.name ==~ /.*\.(cpp|hpp|h)/ || f.name == 'Makefile' })
        files.keySet().sort().each { k ->
            md.update(k.getBytes('UTF-8')); md.update(0 as byte)
            md.update(files[k].bytes);      md.update(0 as byte)
        }
        return md.digest().encodeHex().toString()
    } catch (Exception e) { return "n/a (${e.message})" }
}

// `git describe` of the code directory when it is a git checkout, else n/a.
def gitDescribe() {
    if (!new File(projectDir.toFile(), '.git').exists()) return 'n/a'
    try {
        def p = ['git', '-C', projectDir.toString(), 'describe', '--tags', '--always', '--dirty'].execute()
        p.waitForOrKill(10000)
        def out = p.text.trim()
        return (p.exitValue() == 0 && out) ? out : 'n/a'
    } catch (Exception _e) { return 'n/a' }
}

// ---------------------------------------------------------------------------
// Execution settings (banner and run manifest)
// ---------------------------------------------------------------------------
// The executors and container engines that -profile selects. The workflow
// stops when -profile selects more than one of either.
def profileParts() {
    def profiles = workflow.profile.tokenize(',')*.trim()
    return [ executors: profiles.findAll { p -> p in ['local', 'standard', 'slurm'] },
             engines:   profiles.findAll { p -> p in ['singularity', 'apptainer', 'docker', 'podman'] } ]
}

// 'standard' (no -profile) is the local executor.
def executorName() {
    def e = profileParts().executors
    return (!e || e[0] == 'standard') ? 'local' : e[0]
}

def engineName() {
    def e = profileParts().engines
    return e ? e[0] : 'none'
}

def memUnit(v) {
    return (v instanceof nextflow.util.MemoryUnit) ? v : new nextflow.util.MemoryUnit(v.toString())
}

// The CPUs or memory that a process requests in its first attempt, read from
// the process settings (conf/base.config), and the cap that the executor
// applies: the local limit (conf/local.config) or --max_cpus/--max_memory.
def requestText(String process, String key) {
    def v = nextflow.Global.session.config.process?.get("withName:${process}".toString())?.get(key)
    if (v instanceof Closure) {
        def cl = v.clone()
        cl.delegate = [task: [attempt: 1], params: params]
        cl.resolveStrategy = Closure.DELEGATE_FIRST
        v = cl.call()
    }
    if (v == null) return 'n/a'
    def local = executorName() == 'local'
    def req   = key == 'cpus' ? (v as int) : memUnit(v)
    def cap   = key == 'cpus' ? (local ? params.local_limit_cpus : params.max_cpus)
                              : (local ? params.local_limit_memory : params.max_memory)
    if (cap == null) return "${req}"
    cap = key == 'cpus' ? (cap as int) : memUnit(cap)
    def name = local ? 'local limit' : "max_${key}"
    return req > cap ? "${req} (capped at ${cap}: ${name})" : "${req}"
}

def timeText(t) {
    return executorName() == 'local' ? 'n/a (local)' : "${t}"
}

// Rows shared by the banner and run_manifest.txt. A setting that this run does
// not use shows n/a and the reason.
def stage1Rows() {
    def kmc = params.input_mode == 'kmc'
    return [
        kmer_count_cpus:         kmc ? 'n/a (kmc mode)' : requestText('KMER_COUNT', 'cpus'),
        kmc_to_pack_cpus:        kmc ? requestText('KMC_TO_PACK', 'cpus') : 'n/a (reads mode)',
        kmer_count_memory:       requestText(kmc ? 'KMC_TO_PACK' : 'KMER_COUNT', 'memory'),
        kmer_count_budget_gb:    params.kmer_count_budget_gb == 0 ? '0 (auto)' : "${params.kmer_count_budget_gb} GB",
        kmer_count_read_threads: kmc ? 'n/a (kmc mode)' : params.kmer_count_read_threads,
        kmer_count_time:         timeText(params.kmer_count_time)
    ]
}

def stage2Rows() {
    def off = 'n/a (compress_kbin_packs false)'
    return [
        matrix_merge_cpus:     requestText('MATRIX_MERGE', 'cpus'),
        matrix_merge_memory:   requestText('MATRIX_MERGE', 'memory'),
        matrix_merge_time:     timeText(params.matrix_merge_time),
        compress_packs_cpus:   params.compress_kbin_packs ? requestText('COMPRESS_PACKS', 'cpus') : off,
        compress_packs_memory: params.compress_kbin_packs ? requestText('COMPRESS_PACKS', 'memory') : off
    ]
}

def executionRows() {
    def local = executorName() == 'local'
    def na    = 'n/a (local)'
    return [
        local_limit:           local ? "${params.local_limit_cpus} CPUs / ${memUnit(params.local_limit_memory)}" : 'n/a (slurm)',
        clusterOptions:        local ? na : (params.clusterOptions ?: '(none)'),
        queue_size:            local ? na : params.queue_size,
        max_retries:           local ? na : params.max_retries,
        max_cpus:              params.max_cpus,
        max_memory:            memUnit(params.max_memory),
        max_time:              local ? na : params.max_time,
        singularity_cache_dir: engineName() in ['singularity', 'apptainer'] ? params.singularity_cache_dir
                                                                            : 'n/a (no singularity/apptainer)'
    ]
}

def containerText() {
    def c = workflow.container
    if (!c) return 'none'
    return (c instanceof Map) ? c.values().unique().join(', ') : c.toString()
}

// data_dir is the input only in reads mode; kmc mode reads kmc_dir.
def dataDirText(String data_dir) {
    return params.input_mode == 'kmc' ? 'n/a (kmc mode)' : data_dir
}

// Logs a short summary of pack_stats.tsv: size and coverage ranges, and the
// accessions with the fewest distinct k-mers (often failed or low-coverage samples).
def logPackSummary(f) {
    def lines = f.readLines().findAll { l -> l.trim() }
    if (lines.size() < 2) return
    def keys = lines[0].split('\t') as List
    def rows = lines.drop(1).collect { l -> [keys, l.split('\t', -1) as List].transpose().collectEntries() }
    def nums = { key -> rows.collect { r -> r[key] }.findAll { v -> v && v != 'NA' }.collect { v -> v as long }.sort() }
    def range = { v -> v ? "min ${v[0]}, median ${v[v.size().intdiv(2)]}, max ${v[-1]}" : 'n/a' }
    def lowest = rows.findAll { r -> r.distinct_kmers_after_filter && r.distinct_kmers_after_filter != 'NA' }
                     .sort { r -> (r.distinct_kmers_after_filter as long) }.take(3)
                     .collect { r -> "${r.accession} (${r.distinct_kmers_after_filter})" }.join(', ')
    log.info "Pack statistics (${rows.size()} packs, ${f}):\n" +
             "  distinct k-mers : ${range.call(nums.call('distinct_kmers_after_filter'))}\n" +
             "  k-mer coverage  : ${range.call(nums.call('estimated_kmer_coverage'))} (estimated)\n" +
             "  fewest k-mers   : ${lowest}"
}

// ---------------------------------------------------------------------------
// Run manifest
// ---------------------------------------------------------------------------
// Records the settings a results directory was produced with. Each pack records
// how it was made (input mode, count filter, versions, statistics; see kbin_dump --info),
// but the Stage 2 settings (threshold, encoding, ...) cannot be read back from the
// outputs, so the manifest is the record of how a given results/ directory was made.
def writeManifest(String accessions_file, String data_dir, int n_accessions) {
    def dir = file(params.output_dir)
    dir.mkdirs()
    def mounts = engineName() == 'none' ? 'n/a (no container)' : params.container_mount_dirs.join(', ')
    def manifestRows = (stage1Rows() + stage2Rows() + executionRows() + [container_mount_dirs: mounts])
        .collect { k, v -> "        ${k.padRight(19)} = ${v}" }.join('\n')
    file("${params.output_dir}/run_manifest.txt").text = """\
        # kmer-GWAS run manifest
        # Written at launch; the completion section is added when the run ends.
        pipeline_version    = ${workflow.manifest.version ?: 'n/a'}
        run_name            = ${workflow.runName}
        started             = ${workflow.start}
        status              = running
        nextflow_version    = ${workflow.nextflow.version}
        pipeline_revision   = ${workflow.revision ?: 'n/a'}
        commit_id           = ${workflow.commitId ?: 'n/a'}
        git_describe        = ${gitDescribe()}
        code_checksum       = ${codeChecksum()}
        project_dir         = ${workflow.projectDir}
        launch_dir          = ${workflow.launchDir}
        work_dir            = ${workflow.workDir}
        profile             = ${workflow.profile}
        executor            = ${executorName()}
        container_engine    = ${engineName()}
        container           = ${containerText()}
        command_line        = ${workflow.commandLine}

        kmer_size           = ${kmerSizeText()}
        num_bins            = ${params.num_bins}
        min_kmer_count      = ${minCountText()}
        threshold           = ${params.threshold}
        keep_kmer_counts    = ${params.keep_kmer_counts}
        delimiter           = ${params.delimiter}
        encoding            = ${params.encoding}
        write_core_kmers    = ${params.write_core_kmers}

        accessions_file     = ${accessions_file}
        n_accessions        = ${n_accessions}
        input_mode          = ${params.input_mode}
        kmc_allow_mixed     = ${params.kmc_allow_mixed}
        skip_accession_check = ${params.skip_accession_check}
        kmc_dir             = ${params.input_mode == 'kmc' ? file(params.kmc_dir ?: data_dir).toAbsolutePath().toString() : 'n/a (reads mode)'}
        data_dir            = ${dataDirText(data_dir)}
        output_dir          = ${params.output_dir}

        # Resources and execution. A setting that this run does not use shows n/a.
${manifestRows}

        # Stage 1 packs live in kmer_count_k${params.kmer_size}/. Each pack records
        # its k (${params.kmer_size}) and is rejected by a reader built for any other k.
        # pack_stats.tsv lists the metadata of every pack${params.input_mode == 'kmc' ? ';\n        # kmc_databases.tsv lists the header of every KMC database' : ''}.
        """.stripIndent()
}

// Sets the manifest status and adds the completion section. Only the manifest
// written by this run is changed: a run that stops before writing its own
// manifest leaves an earlier run's manifest as it is.
def writeCompletion() {
    try {
        def f = new File("${params.output_dir}/run_manifest.txt")
        if (!f.isFile()) return
        def text = f.text
        if (!(text =~ ('(?m)^run_name\\s+= ' + java.util.regex.Pattern.quote(workflow.runName) + '$'))) return
        def status = workflow.success ? 'success' : 'failed'
        def err    = workflow.success ? null : (workflow.errorMessage ?: '').readLines().find { l -> l.trim() }?.trim()
        def lines  = [
            '',
            '# Completion',
            "status              = ${status}",
            "completed           = ${workflow.complete}",
            "duration            = ${workflow.duration}",
            "cpu_hours           = ${workflow.stats.computeTimeFmt}",
            "tasks               = ${workflow.stats.succeededCount} succeeded, ${workflow.stats.cachedCount} cached, ${workflow.stats.failedCount} failed",
        ]
        if (!workflow.success) lines << "error               = ${err ?: 'n/a'}"
        f.text = text.replaceFirst(/(?m)^status(\s+)= running$/, "status\$1= ${status}") + lines.join('\n') + '\n'
    } catch (Exception e) { log.warn "Could not add the completion section to run_manifest.txt: ${e.message}" }
}

// ---------------------------------------------------------------------------
// Completion summary
// ---------------------------------------------------------------------------
def completionSummary() {
    writeCompletion()
    // Work-dir cleanup is handled by the native `cleanup` directive in
    // nextflow.config (cleanup = params.cleanup), which is safe for the
    // symlinked publish directories used by MATRIX_MERGE.
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

    if (params.help) {
        helpMessage()
        exit 0
    }

    // Check every parameter against nextflow_schema.json: type, range, allowed
    // values and format. All errors are listed in one message. Rules that involve
    // more than one parameter are checked in the workflow below.
    validateParameters()

    // Backwards compatibility: deprecated flags mapped to current params for
    // continuity. Remove these shims in a future release.
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
    if (params.count != null) {
        log.warn "--count is deprecated and was renamed to --keep_kmer_counts; " +
                 "honouring it as --keep_kmer_counts ${params.keep_kmer_counts}. " +
                 "Please switch, as --count will be removed in a future release."
    }

    workflow.onComplete { completionSummary() }

    // Profiles: at most one executor and one container engine. A container
    // engine needs an executor, because the executor profile sets the resources.
    def executors = profileParts().executors
    def engines   = profileParts().engines
    if (executors.size() > 1)
        exit 1, "ERROR: -profile ${workflow.profile} selects more than one executor (${executors.join(', ')}). Select one: local or slurm."
    if (engines.size() > 1)
        exit 1, "ERROR: -profile ${workflow.profile} selects more than one container engine (${engines.join(', ')}). Select one."
    if (engines && !executors)
        exit 1, "ERROR: -profile ${workflow.profile} has no executor. Add local or slurm, e.g. -profile local,${engines[0]}."

    // The packs record the version in src/version.hpp; it must match the manifest.
    def version_hpp = new File("${params.src_dir}/version.hpp")
    def src_version = null
    if (version_hpp.exists()) {
        def m = version_hpp.text =~ /KMER_GWAS_VERSION\s+"([^"]+)"/
        if (m.find()) src_version = m.group(1)
    }
    if (src_version != workflow.manifest.version)
        exit 1, "ERROR: ${version_hpp} gives version ${src_version}, but manifest.version in " +
                "nextflow.config is ${workflow.manifest.version}. They must match."

    // Resolve k. Reads mode: --kmer_size is required. kmc mode: k is fixed by the
    // databases; nextflow.config reads it from the first accession's .kmc_pre into
    // kmc_detected_kmer_size and adopts it when --kmer_size was not given, so an
    // explicit --kmer_size must match it. A database of a different k is caught by
    // its own conversion job, which aborts on a k mismatch.
    // nextflow_schema.json limits k to an odd number from 15 to 63, also when it
    // comes from the databases. k is compiled in as -DKMER_K, so an invalid value
    // would otherwise fail inside every job. An even k lets a sequence be its own
    // reverse complement, and k must fit two 64-bit words.
    if (params.input_mode == 'kmc') {
        def kdir = params.kmc_dir ?: params.data_dir
        def fmt  = params.kmc_detected_format
        if (params.kmc_detected_kmer_size == null && params.kmer_size == null && fmt != null && fmt != 0x200)
            exit 1, "ERROR: cannot read the k-mer size from this KMC database format " +
                    "(0x${Long.toHexString(fmt as long)}; only KMC2 format 0x200 is read). Give --kmer_size."
        def reason = params.kmc_detect_error ?: "no KMC database header read in ${kdir}"
        if (params.kmc_detected_kmer_size == null && params.kmer_size == null)
            exit 1, "ERROR: could not read the k-mer size from the KMC databases. ${reason}. Expected " +
                    "<kmc_dir>/<accession>.kmc_pre and .kmc_suf for the first accession in ${params.accessions_file}."
        if (params.kmc_detected_kmer_size == null)
            log.warn "Could not read the k-mer size from the KMC databases (${reason}); using --kmer_size " +
                     "${params.kmer_size}. A database built at another k will abort in its conversion job."
        else if ((params.kmer_size as Integer) != (params.kmc_detected_kmer_size as Integer))
            exit 1, "ERROR: --kmer_size ${params.kmer_size} does not match the KMC databases " +
                    "(k=${params.kmc_detected_kmer_size}). Omit --kmer_size, or set it to ${params.kmc_detected_kmer_size}."
    } else if (params.kmer_size == null) {
        exit 1, "ERROR: --kmer_size is required (an ODD number between 15 and 63, e.g. --kmer_size 51)."
    }

    // Validate the output-format combinations up front, so a bad combination
    // fails at launch rather than in every MATRIX_MERGE job. Mirrors the C++ checks.
    if (params.encoding == 'bits' && params.keep_kmer_counts == 'y')
        exit 1, "ERROR: --encoding bits is presence/absence only; it cannot be combined with --keep_kmer_counts y."
    if (params.encoding == 'text' && params.delimiter == 'none' && params.keep_kmer_counts == 'y')
        exit 1, "ERROR: --delimiter none cannot be combined with --keep_kmer_counts y (multi-digit counts would be unparseable). Use --delimiter tab or space."

    // For 'link' publish mode, fail here if output_dir and the work dir are on
    // different filesystems: hard links cannot cross filesystems, so this check
    // avoids a mid-run failure with a clearer error and fix.
    if (params.kmer_count_publish_mode == 'link' || params.matrix_publish_mode == 'link') {
        def outDev = null
        def workDev = null
        try {
            def devOf = { p ->
                // The path itself, or its nearest existing ancestor.
                def path = new File(p.toString()).absoluteFile.toPath()
                def f    = (path.nameCount..0).collect { n -> n == 0 ? path.root : path.root.resolve(path.subpath(0, n)) }
                                                .find { q -> java.nio.file.Files.exists(q) }
                f == null ? null : java.nio.file.Files.getAttribute(f, 'unix:dev')
            }
            outDev  = devOf.call(params.output_dir)
            workDev = devOf.call(workflow.workDir)
        } catch (Exception _e) { /* cannot determine filesystems; let Nextflow surface any link error */ }
        if (outDev != null && workDev != null && outDev != workDev)
            exit 1, "ERROR: 'link' publish mode needs --output_dir and the work dir on the same filesystem, " +
                    "but they are on different ones (output_dir=${params.output_dir}, work=${workflow.workDir}). " +
                    "Set the offending mode (--kmer_count_publish_mode / --matrix_publish_mode) to copy, " +
                    "or point --output_dir at the same filesystem as the work dir."
    }

    // With --kmer_count_publish_mode link and --cleanup false, the uncompressed
    // packs remain hard-linked in the work dir, so --compress_kbin_packs does not
    // reclaim space until the work dir is removed. A warning is issued in this case.
    if (params.compress_kbin_packs && params.kmer_count_publish_mode == 'link' && !params.cleanup)
        log.warn "--compress_kbin_packs is on, but with --kmer_count_publish_mode link and --cleanup false " +
                 "the uncompressed packs remain hard-linked in the work dir, so no space is " +
                 "reclaimed until it is removed. Use --cleanup true (or --kmer_count_publish_mode move)."

    // Resolve to absolute paths — relative inputs break inside Singularity
    // containers and SLURM work directories
    def accessions_file = file(params.accessions_file).toAbsolutePath().toString()
    def data_dir        = file(params.data_dir).toAbsolutePath().toString()
    def kmc_dir         = file(params.kmc_dir ?: params.data_dir).toAbsolutePath().toString()
    // Absolute because Stage 2 and COMPRESS_PACKS receive the pack directory as a
    // path string and run inside their own work directories.
    def output_dir      = file(params.output_dir).toAbsolutePath().toString()

    paramSummary(accessions_file, data_dir)

    def n_accessions = file(accessions_file).readLines().findAll { l -> l.trim() }.size()
    writeManifest(accessions_file, data_dir, n_accessions)

    // -- Stage 1: one pack per accession, counted from reads or converted from a
    //    KMC database. Both processes emit the same (accession, pack) tuples. --
    ch_accessions = channel
        .fromPath(params.accessions_file)
        .splitText()
        .map { l -> l.trim() }
        .filter { l -> l }                  // skip blank lines

    // Build the bits_to_text converter into results/bin/ so the output
    // directory ships with the tool needed to read its bit-packed matrices.
    BUILD_TOOLS()

    if (params.input_mode == 'kmc') {
        // Check every database header first; the conversions start only after
        // the check has passed.
        KMC_HEADERS(file(accessions_file), kmc_dir)
        KMC_HEADERS.out.notes
            .splitText().map { l -> l.trim() }.filter { l -> l }
            .subscribe { l -> l.startsWith('WARNING:') ? log.warn(l - 'WARNING: ') : log.info(l - 'NOTE: ') }
        ch_checked = ch_accessions.combine(KMC_HEADERS.out.table).map { acc, _table -> acc }
        KMC_TO_PACK(ch_checked, params.num_bins, kmc_dir)
        ch_packs = KMC_TO_PACK.out.accession_pack
        ch_stats = KMC_TO_PACK.out.stats
    } else {
        KMER_COUNT(ch_accessions, params.num_bins, data_dir)
        ch_packs = KMER_COUNT.out.accession_pack
        ch_stats = KMER_COUNT.out.stats
    }

    // One table with the metadata of every pack, written when Stage 1 is done,
    // with the rows in the order of the accessions file.
    def acc_order = [:]
    file(accessions_file).readLines()*.trim().findAll { l -> l }.eachWithIndex { a, i -> acc_order[a] = i }
    ch_stats
        .collectFile(name: 'pack_stats.tsv', storeDir: output_dir, keepHeader: true,
                     sort: { f -> acc_order.get(f.name - ~/\.pack_stats\.tsv$/, Integer.MAX_VALUE) })
        .subscribe { f -> logPackSummary(f) }

    // Fan-in: wait for all accessions, then fire once per bin. Each finished
    // accession emits num_bins sentinel tuples (bin_idx, kmer_count_root).
    // groupTuple() collects one entry per accession for every bin_idx; the
    // channel only closes (and each bin fires) once all accessions have
    // completed Stage 1.
    // Must match the Stage 1 publishDir exactly (KMER_COUNT and KMC_TO_PACK use
    // the same one) — k is part of the path so that a results directory from a
    // previous run at a different k cannot be picked up silently.
    def kmer_count_root = "${output_dir}/kmer_count_k${params.kmer_size}"
    def n_bins          = params.num_bins as Integer

    ch_bin_signals = ch_packs
        .flatMap { _accession, _pack ->
            (0..<n_bins).collect { bin_idx ->
                tuple(bin_idx, kmer_count_root)
            }
        }

    // One directory per bin, passed as a VALUE (absolute path string), not a
    // Nextflow `path` input. Two reasons: (1) passing the packs individually would
    // cost one staged symlink per accession per bin task (millions of inodes for a
    // large cohort at queueSize 200); (2) the pack directory is the PUBLISHED
    // output_dir, which is mutated after a run — re-publish timestamps, and
    // --compress_kbin_packs rewrites each .kbin to .kbin.gz in place. Hashing it as
    // a `path` input made -resume always re-run MATRIX_MERGE (and, once packs were
    // gzipped, fail). As a value the cache key is the stable path string; correct
    // invalidation still comes from k (in the path), the accessions file, and the
    // Stage-2 params baked into the script (see the signature line in the module).
    ch_bins_ready = ch_bin_signals
        .groupTuple()                           // group by bin_idx; size == num_accessions
        .map { bin_idx, roots -> tuple(bin_idx, roots[0]) }

    // -- Stage 2: matrix merge, one job per bin --
    MATRIX_MERGE(ch_bins_ready, file(accessions_file))

    // -- Optional: compress the published packs, one job per accession, only AFTER
    //    every bin's MATRIX_MERGE has finished (gated on the collected outputs), so
    //    a pack is never compressed while it might still be read. See the module.
    if (params.compress_kbin_packs) {
        ch_gate = MATRIX_MERGE.out.matrix_file.collect().map { true }
        ch_accs = ch_packs.map { acc, _pack -> acc }
        COMPRESS_PACKS(ch_accs.combine(ch_gate), kmer_count_root)
    }
}
