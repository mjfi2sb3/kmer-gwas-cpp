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
    # Load gcc if not available or below minimum version (9)
    if ! command -v g++ &>/dev/null || [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        module load gcc/12.2.0 2>/dev/null || true
    fi
    if ! command -v g++ &>/dev/null; then
        echo "WARNING: g++ not available after module load attempt" >&2
    elif [[ \$(g++ -dumpversion | cut -d. -f1) -lt 9 ]]; then
        echo "WARNING: g++ \$(g++ -dumpversion) < 9 -- compilation may fail" >&2
    fi

    g++ -std=c++17 -O3 -march=native -pthread -o kmer_count \
        ${projectDir}/src/kmer_count_v3.cpp \
        ${projectDir}/src/mmap_io.cpp

    mkdir -p data

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

    # Stage R1 — decompress if needed, otherwise symlink
    if [[ "\$R1" == *.gz ]]; then
        gunzip -c "\$R1" > data/${accession}_1.fq
    else
        ln -sf "\$R1" data/${accession}_1.fq
    fi

    # Stage R2
    if [[ "\$R2" == *.gz ]]; then
        gunzip -c "\$R2" > data/${accession}_2.fq
    else
        ln -sf "\$R2" data/${accession}_2.fq
    fi

    ./kmer_count ${accession} ${num_bins} ./
    """
}
