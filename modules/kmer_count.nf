// Stage 1 — count one accession's k-mers into a single self-indexed pack file.
//
// Emits `<accession>.kbin`, replacing the previous per-bin directory tree plus
// tar archive. Inode cost per task went from 1 + num_bins dirs + 2 x num_bins
// files (measured 4,502 peak at 1500 bins) to one file.
//
// publishDir uses mode 'move', not 'copy': the pack is large (~5 GB per rice
// accession) and copying would keep a second copy in the work directory until
// the end of the run — 63 TB of duplication at 12,600 accessions. Moving is
// safe here because MATRIX_MERGE consumes the published DIRECTORY, not the
// individual files through a channel (see main.nf).
process KMER_COUNT {
    tag "${accession}"

    // k is in the directory name so bin files from different k-mer lengths can
    // never be mixed: they store fixed-width k-mers and are unreadable at any
    // other k. MATRIX_MERGE reads this same path (see main.nf kmer_count_root).
    publishDir "${params.output_dir}/kmer_count_k${params.kmer_size}", mode: 'move', overwrite: true

    input:
        val accession
        val num_bins
        val data_dir

    output:
        tuple val(accession), path("${accession}.kbin"), emit: accession_pack

    script:
    // Give the counter most of the task's memory as its accumulation budget,
    // leaving headroom for the read buffers and the runtime itself.
    // task.memory is null under profiles that set no per-process memory
    // directive (e.g. `standard`), so fall back to a modest default.
    def budget_gb = task.memory ? Math.max(1, (int)(task.memory.toGiga() * 0.7)) : 8
    """
    # -DKMER_K makes the k-mer length a compile-time constant, so it is
    # user-configurable at no runtime cost (the binary is built per job anyway,
    # to get -march=native for this node's CPU).
    g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=${params.kmer_size} -o kmer_count \
        ${params.src_dir}/kmer_count_v3.cpp -lz

    # Locate R1 — try common extensions in order
    R1=""
    for ext in _1.fq _1.fastq _1.fq.gz _1.fastq.gz; do
        candidate="${data_dir}/${accession}\${ext}"
        if [ -f "\${candidate}" ]; then R1="\${candidate}"; break; fi
    done
    [ -z "\$R1" ] && { echo "ERROR: no R1 file found for ${accession} in ${data_dir}" >&2; exit 1; }

    # Locate R2
    R2=""
    for ext in _2.fq _2.fastq _2.fq.gz _2.fastq.gz; do
        candidate="${data_dir}/${accession}\${ext}"
        if [ -f "\${candidate}" ]; then R2="\${candidate}"; break; fi
    done
    [ -z "\$R2" ] && { echo "ERROR: no R2 file found for ${accession} in ${data_dir}" >&2; exit 1; }

    ./kmer_count ${accession} ${num_bins} ./ "\$R1" "\$R2" ${budget_gb} ${params.min_kmer_count}
    """
}
