FROM gcc:12
RUN apt-get update && apt-get install -y zlib1g-dev pigz cmake && rm -rf /var/lib/apt/lists/*

# zlib-ng — a faster gzip inflate (~3x zlib on FASTQ). Built from source because
# Debian has no zlib-ng package. WITH_GZFILEOP exposes the zng_gz* file ops the
# counter uses; native symbols (ZLIB_COMPAT=OFF) leave stock zlib in place as the
# portable fallback. kmer_count probes for it and links whichever is present.
RUN git clone --depth 1 --branch 2.3.3 https://github.com/zlib-ng/zlib-ng.git /tmp/zlib-ng \
    && cmake -S /tmp/zlib-ng -B /tmp/zlib-ng/build \
        -DCMAKE_BUILD_TYPE=Release -DWITH_GZFILEOP=ON -DZLIB_COMPAT=OFF -DBUILD_SHARED_LIBS=ON \
    && cmake --build /tmp/zlib-ng/build -j"$(nproc)" \
    && cmake --install /tmp/zlib-ng/build \
    && ldconfig \
    && rm -rf /tmp/zlib-ng

COPY src/ /opt/kmer-gwas/src/
