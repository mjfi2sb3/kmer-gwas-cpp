#ifndef VERSION_HPP
#define VERSION_HPP

// Pipeline version, recorded in the metadata of every pack. It must equal
// manifest.version in nextflow.config; main.nf stops if they differ.
//
// It is kept here, not passed on the compile line, so that the Nextflow task
// scripts do not change between releases and -resume keeps its cached tasks.
#define KMER_GWAS_VERSION "3.7.5"

#endif
