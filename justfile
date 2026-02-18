set shell := ["bash", "-cu"]

# Delete all cmake-build-
clear:
    rm -rf cmake-build-*

# Baseline (no AVX) build
baseline:
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target protal -- -j$(nproc)

# AVX2 build
avx2:
    cmake -S . -B cmake-build-avx2 -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-avx2 --target protal_avx2 -- -j$(nproc)

# Baseline simulate_metagenomes build
simulate:
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target simulate_metagenomes -- -j$(nproc)

# Static baseline build (protal_static target)
static:
    cmake -S . -B cmake-build-static -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-static --target protal_static -- -j$(nproc)

# Build all targets, baseline avx2 static and simulate
build-all: clear baseline avx2 static simulate
