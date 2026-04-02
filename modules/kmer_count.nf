process KMER_COUNT {
    tag "${accession}"
    publishDir "${params.output_dir}/kmer_count", mode: params.large_cohort ? 'move' : 'copy', overwrite: true

    input:
        val accession
        val num_bins
        val data_dir

    output:
        tuple val(accession), path(params.large_cohort ? "${accession}.tar" : "${accession}"), emit: accession_dir

    script:
    def pack = params.large_cohort ? "tar -cf ${accession}.tar ${accession}/ && rm -rf ${accession}/" : ""
    """
    g++ -std=c++17 -O3 -march=native -pthread -o kmer_count \
        /opt/kmer-gwas/src/kmer_count_v3.cpp \
        /opt/kmer-gwas/src/mmap_io.cpp -lz

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

    ./kmer_count ${accession} ${num_bins} ./ "\$R1" "\$R2"

    ${pack}
    """
}
