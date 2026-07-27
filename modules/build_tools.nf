// Compile the bits_to_text converter and publish it alongside the results, so
// a run's output directory is self-contained: the matrices plus the tool that
// converts the bit-packed form to and from text.
//
// Unlike the per-job kmer_count/matrix_merge binaries (built with -march=native
// for the compute node), this is built WITHOUT -march=native: it is a small
// I/O-bound utility the user may run anywhere — a login node, a laptop — so
// portability across CPUs matters more than the negligible speed-up.
process BUILD_TOOLS {
    publishDir "${params.output_dir}/bin", mode: 'copy', overwrite: true

    output:
        path 'bits_to_text', emit: bits_to_text

    script:
    """
    g++ -std=c++17 -O2 -o bits_to_text ${params.src_dir}/bits_to_text.cpp -lz
    """
}
