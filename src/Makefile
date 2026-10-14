# Local build for the pipeline binaries: kmer_count, matrix_merge, bits_to_text,
# kbin_dump and kmc_to_pack. kmc_to_pack links the KMC C++ API in kmc_api/
# (GNU GPL 3), so that binary is GPL.
#
# In production the Nextflow modules compile these at job start (to get
# -march=native for the actual compute node), so this Makefile exists for local
# development and testing. It mirrors the flags the modules use.
#
#   make                 # build all at the default k
#   make KMER_K=31       # build for 31-mers (k <= 32 uses the compact 10-byte record)
#   make test            # build, then run a tiny end-to-end sanity check of
#                        # kmer_count + matrix_merge (kmc_to_pack is tested by
#                        # ../tools/test_kmc_to_pack.sh, which needs KMC)
#   make clean

KMER_K  ?= 51
CXX     ?= g++
CXXFLAGS = -std=c++17 -O3 -march=native -pthread -DKMER_K=$(KMER_K) -Wall -Wextra
BIN     ?= .

HEADERS = kmer_key.hpp pack_io.hpp bin_store.hpp thread_pool.hpp mem_limit.hpp version.hpp

all: $(BIN)/kmer_count $(BIN)/matrix_merge $(BIN)/bits_to_text $(BIN)/kbin_dump $(BIN)/kmc_to_pack

$(BIN)/kmer_count: kmer_count_v3.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ kmer_count_v3.cpp -lz

$(BIN)/matrix_merge: matrix_merge.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ matrix_merge.cpp

# Standalone; does not depend on k or the key format (a pure text transform).
$(BIN)/bits_to_text: bits_to_text.cpp
	$(CXX) $(CXXFLAGS) -o $@ bits_to_text.cpp -lz

# Standalone; reads the pack's k from its footer, so one binary reads any pack.
$(BIN)/kbin_dump: kbin_dump.cpp
	$(CXX) $(CXXFLAGS) -o $@ kbin_dump.cpp -lz

# Converts a KMC database into a .kbin pack, reading it through the vendored KMC
# C++ API (src/kmc_api/, GNU GPL 3 — so this binary is GPL). k baked in via -DKMER_K.
KMC_API = kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp
$(BIN)/kmc_to_pack: kmc_to_pack.cpp $(HEADERS) $(KMC_API)
	$(CXX) $(CXXFLAGS) -o $@ kmc_to_pack.cpp $(KMC_API)

# Round-trip one synthetic accession through both stages.
test: all
	@set -e; tmp=$$(mktemp -d); \
	python3 -c "import random;random.seed(0);\
g=''.join(random.choice('ACGT') for _ in range(20000));\
fh=open('$$tmp/T_1.fq','w');\
[fh.write('@r%d\n%s\n+\n%s\n'%(i,g[(i*7)%(len(g)-150):(i*7)%(len(g)-150)+150],'I'*150)) for i in range(2000)];\
fh.close();import shutil;shutil.copy('$$tmp/T_1.fq','$$tmp/T_2.fq')"; \
	echo T > $$tmp/acc.txt; \
	(cd $$tmp && $(abspath $(BIN))/kmer_count T 4 ./ T_1.fq T_2.fq 1 2 > /dev/null); \
	(cd $$tmp && $(abspath $(BIN))/matrix_merge --input ./ --accessions acc.txt --index 0 \
	    --threshold 0 --delimiter tab --count n --core n --bins 4 --threads 1 > /dev/null); \
	rows=$$(cat $$tmp/matrix_*/0_matrix.tsv | wc -l); \
	width=$$(head -1 $$tmp/matrix_*/0_matrix.tsv | cut -f1 | tr -d '\n' | wc -c); \
	echo "  k=$(KMER_K): $$rows rows, k-mer column width $$width (expect $(KMER_K))"; \
	test "$$width" = "$(KMER_K)" || { echo "  FAIL: wrong k-mer width"; exit 1; }; \
	rm -rf $$tmp; echo "  OK"

clean:
	rm -f $(BIN)/kmer_count $(BIN)/matrix_merge $(BIN)/bits_to_text $(BIN)/kbin_dump $(BIN)/kmc_to_pack

.PHONY: all test clean
