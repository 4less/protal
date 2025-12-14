set shell := ["bash", "-cu"]

# Baseline (no AVX) build
baseline:
    rm -rf cmake-build-baseline
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target protal -- -j$(nproc)

# AVX2 build
avx2:
    rm -rf cmake-build-avx2
    cmake -S . -B cmake-build-avx2 -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-avx2 --target protal_avx2 -- -j$(nproc)

# Baseline simulate_metagenomes build
simulate:
    rm -rf cmake-build-baseline
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target simulate_metagenomes -- -j$(nproc)

# Static baseline build using musl (no AVX)
musl-static:
    MUSL_TC="${MUSL_TOOLCHAIN:-$HOME/musl-toolchain/x86_64-linux-musl-native/bin}"; \
    [ -x "${MUSL_TC}/x86_64-linux-musl-g++" ] || { echo "musl toolchain not found at ${MUSL_TC}; set MUSL_TOOLCHAIN or install the musl.cc toolchain there." >&2; exit 1; }; \
    rm -rf cmake-build-musl; \
    CC="${MUSL_TC}/x86_64-linux-musl-gcc" CXX="${MUSL_TC}/x86_64-linux-musl-g++" \
    CFLAGS="-O3 -pthread -march=x86-64 -mtune=generic -mno-avx -mno-avx2 -fno-tree-vectorize -static -static-libgcc -Wno-error=unknown-pragmas" \
    CXXFLAGS="-O3 -pthread -march=x86-64 -mtune=generic -mno-avx -mno-avx2 -fno-tree-vectorize -static -static-libgcc -static-libstdc++ -Wno-error=unknown-pragmas" \
    LDFLAGS="-static -static-libgcc -static-libstdc++ -L/opt/musl-openmp/x86_64-linux-musl/lib" \
    CMAKE_PREFIX_PATH="/opt/musl-openmp/x86_64-linux-musl" \
    PKG_CONFIG_PATH="/opt/musl-openmp/x86_64-linux-musl/lib/pkgconfig" \
    cmake -S . -B cmake-build-musl \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SHARED_LIBS=OFF \
      -DZLIB_LIBRARY=/opt/musl-openmp/x86_64-linux-musl/lib/libz.a \
      -DZLIB_INCLUDE_DIR=/opt/musl-openmp/x86_64-linux-musl/include; \
    cmake --build cmake-build-musl --target protal -- -j$(nproc)


# Static baseline build using musl (no AVX)
musl-static-avx2:
    MUSL_TC="${MUSL_TOOLCHAIN:-$HOME/musl-toolchain/x86_64-linux-musl-native/bin}"; \
    [ -x "${MUSL_TC}/x86_64-linux-musl-g++" ] || { echo "musl toolchain not found at ${MUSL_TC}; set MUSL_TOOLCHAIN or install the musl.cc toolchain there." >&2; exit 1; }; \
    rm -rf cmake-build-musl; \
    CC="${MUSL_TC}/x86_64-linux-musl-gcc" CXX="${MUSL_TC}/x86_64-linux-musl-g++" \
    CFLAGS="-O3 -pthread -march=x86-64 -mtune=generic -static -static-libgcc -Wno-error=unknown-pragmas" \
    CXXFLAGS="-O3 -pthread -march=x86-64 -mtune=generic -static -static-libgcc -static-libstdc++ -Wno-error=unknown-pragmas" \
    LDFLAGS="-static -static-libgcc -static-libstdc++ -L$HOME/musl-toolchain/x86_64-linux-musl-native/lib" \
    CMAKE_PREFIX_PATH="$HOME/musl-toolchain/x86_64-linux-musl-native/" \
    PKG_CONFIG_PATH="$HOME/musl-toolchain/x86_64-linux-musl-native/lib/pkgconfig" \
    cmake -S . -B cmake-build-musl-avx2 \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SHARED_LIBS=OFF \
      -DZLIB_LIBRARY="$HOME/musl-toolchain/x86_64-linux-musl-native/lib/libz.a" \
      -DZLIB_INCLUDE_DIR="$HOME/musl-toolchain/x86_64-linux-musl-native/include"; \
    cmake --build cmake-build-musl-avx2 --target protal_avx2 -- -j$(nproc)

# Build both
build-all: baseline avx2
