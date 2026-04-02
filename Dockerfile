FROM gcc:12
RUN apt-get update && apt-get install -y zlib1g-dev pigz && rm -rf /var/lib/apt/lists/*
COPY src/ /opt/kmer-gwas/src/
