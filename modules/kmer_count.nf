process KMER_COUNT {
    tag "${accession}"
    publishDir "${params.output_dir}/kmer_count", mode: 'copy', overwrite: true

    input:
        val accession
        val num_bins

    output:
        tuple val(accession), path("${accession}"), emit: accession_dir

    script:
    """
    # Ensure module system is available (not sourced in non-login shells)
    source /etc/profile.d/modules.sh 2>/dev/null || true

    # Load gcc if not available or below minimum version (9)
    if ! command -v g++ &>/dev/null || [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        module load gcc/12.2.0 2>/dev/null || true
    fi
    if ! command -v g++ &>/dev/null; then
        echo "WARNING: g++ not available after module load attempt" >&2
    elif [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        echo "WARNING: g++ \$(g++ -dumpversion) < 9 -- compilation may fail" >&2
    fi

    module load zlib/1.2.13/gnu-12.2.0 2>/dev/null || true

    g++ -std=c++17 -O3 -march=native -pthread -o kmer_count \
        ${projectDir}/src/kmer_count_v3.cpp \
        ${projectDir}/src/mmap_io.cpp -lz

    # Locate R1 — try common extensions in order
    R1=""
    for ext in _1.fq _1.fastq _1.fq.gz _1.fastq.gz; do
        candidate="${params.data_dir}/${accession}\${ext}"
        if [ -f "\${candidate}" ]; then R1="\${candidate}"; break; fi
    done
    [ -z "\$R1" ] && { echo "ERROR: no R1 file found for ${accession} in ${params.data_dir}" >&2; exit 1; }

    # Locate R2
    R2=""
    for ext in _2.fq _2.fastq _2.fq.gz _2.fastq.gz; do
        candidate="${params.data_dir}/${accession}\${ext}"
        if [ -f "\${candidate}" ]; then R2="\${candidate}"; break; fi
    done
    [ -z "\$R2" ] && { echo "ERROR: no R2 file found for ${accession} in ${params.data_dir}" >&2; exit 1; }

    ./kmer_count ${accession} ${num_bins} ./ "\$R1" "\$R2"
    """
}
