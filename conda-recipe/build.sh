#!/bin/bash
echo $PWD
ls $PWD

mkdir build
cd build


echo "Step 1"
echo $PWD


cmake \
	-DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_MAKE_PROGRAM=make \
	-DCMAKE_C_COMPILER=$GCC \
	-DCMAKE_CXX_COMPILER=$GXX \
	-G "CodeBlocks - Unix Makefiles" \
	-S ./../ \
	-B cmake-build-release \
        -DCMAKE_VERBOSE_MAKEFILE=ON

echo "Step 2"
echo $PWD
echo "cmake --build cmake-build-release --target protal -- -j 6"

cmake --build cmake-build-release --target protal -- -j 6
cmake --build cmake-build-release --target protal_avx2 -- -j 6

cp cmake-build-release/protal ${PREFIX}/bin/protal_plain
cp cmake-build-release/protal_avx2 ${PREFIX}/bin/protal_avx2
cp ./../protal_launcher ${PREFIX}/bin/protal
