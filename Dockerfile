FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libgmp-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Clone and build OpenFHE v1.1.2
RUN git clone --branch v1.1.2 --depth 1 https://github.com/openfhe-development/openfhe-development.git && \
    cd openfhe-development && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf openfhe-development

# Set library path
ENV LD_LIBRARY_PATH=/usr/local/lib

# Copy project files
COPY data /app/data
COPY src /app/src
COPY benchmark.sh /app/benchmark.sh
COPY CMakeLists.txt /app/CMakeLists.txt
COPY LICENSE /app/LICENSE
COPY PreLoad.cmake /app/PreLoad.cmake
COPY README.md /app/README.md

# Build the project
RUN mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc)

# Set default command
CMD ["/bin/bash"]
