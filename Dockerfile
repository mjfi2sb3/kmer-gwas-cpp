FROM gcc:12
RUN apt-get update && apt-get install -y zlib1g-dev && rm -rf /var/lib/apt/lists/*
