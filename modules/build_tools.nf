// Compile the helper tools and publish them alongside the results, so a run's
// output directory is self-contained: the matrices and packs plus the tools that
// read them.
//
//   bits_to_text — convert the bit-packed matrix to and from text.
//   kbin_dump    — inspect a Stage 1 .kbin pack or export its k-mers as text.
//
// Unlike the per-job kmer_count/matrix_merge binaries (built with -march=native
// for the compute node), these are built WITHOUT -march=native: they are small
// I/O-bound utilities the user may run anywhere — a login node, a laptop — so
// portability across CPUs matters more than the negligible speed-up. kbin_dump
// reads a pack's k from its footer, so one binary reads a pack of any k.
process BUILD_TOOLS {
    publishDir "${params.output_dir}/bin", mode: 'copy', overwrite: true

    output:
        path 'bits_to_text', emit: bits_to_text
        path 'kbin_dump',    emit: kbin_dump

    script:
    """
    g++ -std=c++17 -O2 -o bits_to_text ${params.src_dir}/bits_to_text.cpp -lz
    g++ -std=c++17 -O2 -o kbin_dump    ${params.src_dir}/kbin_dump.cpp
    """
}
