// Stage 1 — count one accession's k-mers into a single self-indexed pack file.
//
// Emits `<accession>.kbin`, replacing the previous per-bin directory tree plus
// tar archive. Inode cost per task went from 1 + num_bins dirs + 2 x num_bins
// files (measured 4,502 peak at 1500 bins) to one file.
//
// The pack is large (~5 GB per rice accession), so the publish mode is a real
// trade-off, exposed as --publish_mode (see nextflow.config):
//   'link' (default) — a hard link: no second copy of the data, and the pack
//     stays in the work dir so -resume can skip finished accessions. Needs
//     output_dir and the work dir on one filesystem (main.nf checks this).
//   'copy' — a second copy in output_dir; -resume works on any layout but keeps
//     ~2x storage until cleanup.
//   'move' — no second copy, but the pack leaves the work dir so -resume cannot
//     skip finished accessions.
// MATRIX_MERGE consumes the published DIRECTORY, not the packs through a channel
// (see main.nf), so any of these modes works for the dataflow.
process KMER_COUNT {
    tag "${accession}"

    // k is in the directory name so bin files from different k-mer lengths can
    // never be mixed: they store fixed-width k-mers and are unreadable at any
    // other k. MATRIX_MERGE reads this same path (see main.nf kmer_count_root).
    publishDir "${params.output_dir}/kmer_count_k${params.kmer_size}", mode: params.publish_mode, overwrite: true

    input:
        val accession
        val num_bins
        val data_dir

    output:
        tuple val(accession), path("${accession}.kbin"), emit: accession_pack

    script:
    // The accumulation budget is auto-sized INSIDE kmer_count from the memory
    // actually enforced on the job (cgroup / SLURM env / MemTotal), so it tracks
    // the real node — a shared allocation, an --mem=0 exclusive node, or a bare
    // workstation — without depending on task.memory (which is just the SLURM
    // request and is wrong under --mem=0). Passing 0 selects that AUTO sizing;
    // --kmer_count_budget_gb sets an explicit budget, capped at the enforced
    // limit. See src/mem_limit.hpp.
    """
    # Prefer zlib-ng for gzip inflate (~3x faster than stock zlib on FASTQ) when
    # it is available (it is in the container image); fall back to stock zlib so
    # the non-container profiles still build on clusters without it. Probe by
    # actually compiling+linking a tiny program against zlib-ng's gz ops.
    ZFLAGS="-lz"
    if printf '#define WITH_GZFILEOP\\n#include <zlib-ng.h>\\nint main(){ gzFile f=zng_gzopen("","rb"); (void)f; return 0; }\\n' \
         | g++ -x c++ - -lz-ng -o /dev/null 2>/dev/null; then
        ZFLAGS="-DHAVE_ZLIBNG -lz-ng"
    fi

    # -DKMER_K makes the k-mer length a compile-time constant, so it is
    # user-configurable at no runtime cost (the binary is built per job anyway,
    # to get -march=native for this node's CPU).
    g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=${params.kmer_size} \$ZFLAGS -o kmer_count \
        ${params.src_dir}/kmer_count_v3.cpp

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

    ./kmer_count ${accession} ${num_bins} ./ "\$R1" "\$R2" ${params.kmer_count_budget_gb} ${params.min_kmer_count} ${params.kmer_count_read_threads}
    """
}
