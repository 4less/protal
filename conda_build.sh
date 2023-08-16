if [ -z ${CONDA_BUILD+x} ]; then
    source /usr/users/QIB_fr017/fritsche/miniconda3/envs/conda-install-test/conda-bld/protal_1692176526483/work/build_env_setup.sh
fi
#!/bin/bash

echo "------------------------------------------------------------------"
echo ""
echo ""
echo ""
echo $PWD
ls -1

which g++
which cc
which ninja


rm -r cmake-build-release
mkdir cmake-build-release

# CMAKE_PATH=cmake
CC_PATH=$BUILD_PREFIX/bin/cc
GPP_PATH=$BUILD_PREFIX/bin/x86_64-conda_cos7-linux-gnu-g++
# NINJA_PATH=$BUILD_PREFIX/bin/ninja

export CMAKE_C_COMPILER=$BUILD_PREFIX/bin/gcc
export CMAKE_CXX_COMPILER=$BUILD_PREFIX/bin/g++
export CMAKE_INCLUDE_PATH=$BUILD_PREFIX/include
export CMAKE_LIBRARY_PATH=$BUILD_PREFIX/lib

CMAKE_PATH=/usr/users/QIB_fr017/fritsche/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/231.9225.21/bin/cmake/linux/x64/bin/cmake
# CC_PATH=/usr/bin/cc
# GPP_PATH=/usr/bin/g++
NINJA_PATH=/usr/users/QIB_fr017/fritsche/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/231.9225.21/bin/ninja/linux/x64/ninja

LD_LIBRARY_PATH=$BUILD_PREFIX/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
$CMAKE_PATH \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_C_COMPILER=$CC_PATH \
	-DCMAKE_CXX_COMPILER=$GPP_PATH \
	-DBUILD_SHARED_LIBS=OFF \
	-S ./ \
	-B cmake-build-release \
	-G Ninja \
	-DCMAKE_MAKE_PROGRAM=$NINJA_PATH




ls -1 cmake-build-release
# cd cmake-build-release
$NINJA_PATH -C cmake-build-release
# $BUILD_PREFIX/bin/make
# cp cmake-build-release/protal $PREFIX/bin/protal


echo ""
echo ""
echo ""
echo "------------------------------------------------------------------"