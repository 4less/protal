project(protal_lib)
set(CMAKE_CXX_STANDARD 20)

# -lcPMML
# set(CMAKE_EXE_LINKER_FLAGS " -lwfa -flto") #these flags are not needed because set somewhere else
# set(CMAKE_CXX_FLAGS_DEBUG " -pthread -ggdb -g -std=c++2a -fopenmp -Wall -mavx2 -mavx -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx -march=native") #-march=native
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx -march=native") #-march=native

# CHECK DEBUG
if(DEFINED ENV{GCC} AND NOT "$ENV{GCC}" STREQUAL "")
    message(STATUS "Using GCC-prefixed binutils from Conda: $ENV{GCC}")
    message(STATUS "Using GCC Prefix: ${GCC_PREFIX}")
    message(STATUS "PREFIX PATH: ${PREFIX}")
    SET(CMAKE_AR  "$ENV{GCC}-ar")
    SET(CMAKE_RANLIB  "$ENV{GCC}-ranlib")
else()
    message(STATUS "Using system default ar/ranlib")
    set(CMAKE_AR "ar")
    set(CMAKE_RANLIB "ranlib")
endif()


SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_C_ARCHIVE_FINISH   true)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")

message(${CMAKE_CXX_FLAGS_DEBUG})
message("CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")



#set(CMAKE_PREFIX_PATH "/usr/users/QIB_fr017/fritsche/miniconda3/envs/gcconda/" ${CMAKE_PREFIX_PATH})

# This creates the target and must be specified before using target_include etc
add_library(protal_lib
        RunProtal.h
        Options.h
        Build.h
        Classify.h
        Utils/Benchmark.h
        Utils/Utils.h
        Utils/ProgressBar.h
        Utils/xxhash64.h
        Utils/BinaryClassifierEvaluator.h
        Utils/LineSplitter.h
        SequenceUtils/SeqReader.h
        SequenceUtils/KmerIterator.h
        SequenceUtils/KmerUtils.h
        SequenceUtils/Minimizer.h
        Hash/compact_map.h
        Hash/compact_map.cpp
        Hash/compact_map_utils.h
        Hash/Seedmap.h
        Hash/SeedmapFlex.h
        Hash/KmerPutter.h
        Hash/KmerLookup.h
        Profiling/Profiler.h
        Profiling/Strain.h
        Profiling/Coverage.h
        Profiling/Coverage.cpp
        Profiling/Coverage.h
        Profiling/ProfilerDefinitions.h
        SequenceUtils/GenomeLoader.h
        IO/AlignmentOutputHandler.h
        Alignment/WFA2Wrapper.h
        Alignment/FastAlignment.h
        Alignment/InternalReadAlignment.h
        Alignment/ScoreAlignments.h
        Alignment/AlignmentUtils.h
        Core/CoreBenchmark.h
        IO/FastxReader.h
        IO/FastxReader.cpp
        IO/SamHandler.h
        IO/BufferedOutput.h
        IO/BufferedOutput.cpp
        Core/SeedingStrategy.h
        Core/ChainingStrategy.h
        Core/AnchorFinder.h
        Core/ChainAnchorFinder.h
        Core/AlignmentStrategy.h
        VarkitInterface/VarkitInterface.h
        Taxonomy/Taxonomy.h
        Taxonomy/Taxonomy.cpp
        Statistics.h
        Constants.h
        SNPs/SNP.h
        SNPs/SNPUtils.h
        Test.h
        SequenceUtils/SequenceRange.h
        SequenceUtils/SequenceRange.cpp
        SequenceUtils/SequenceRangeHandler.h
        SequenceUtils/SequenceRangeHandler.cpp
        SequenceUtils/Variant.h
        SequenceUtils/VariantHandler.h
        Strains/Distances.h
        IO/TaxonStatisticsOutput.h
        Utils/Matrix.h
        Strains/MSA.h
        Profiling/TreeFilter.h
        Alignment/WFA2Wrapper2.h
        )

target_include_directories(protal_lib
        PUBLIC
        ../lib
        IO
        Hash
        Core
        VarkitInterface
        Alignment
        Utils
        Profiling
        SequenceUtils
        Taxonomy
        SNPs)


set_target_properties(protal_lib PROPERTIES LINKER_LANGUAGE CXX)

target_link_directories(protal_lib
        PUBLIC
        ../IO
        ../lib
        ../Hash
        ../Core
        ../VarkitInterface
        ../Alignment
        ../Utils
        ../Profiling
        # ../Profiler
        ../SequenceUtils
        ../Taxonomy
        ../SNPs
        ../cPMML
)

# target_include_directories(protal_lib
#         PUBLIC
#         ../../lib
#         ../IO
#         ../Hash
#         ../Core
#         ../VarkitInterface
#         ../Alignment
#         ../Utils
#         ../Profiling
#         ../SequenceUtils
#         ../Taxonomy
#         ../SNPs)


#set(ZLIB_USE_STATIC_LIBS "ON")
find_package(ZLIB REQUIRED)

target_include_directories(protal_lib PUBLIC ../lib/robin)
target_include_directories(protal_lib PUBLIC ../lib/tsl)

# Add-Subdirectory with ${CMAKE_CURRENT_SOURCE_DIR}/


add_subdirectory(${PROJECT_SOURCE_DIR}/../lib/wfa2-lib wfa)
add_subdirectory(${PROJECT_SOURCE_DIR}/../lib/gzstream gzstream)
add_subdirectory(${PROJECT_SOURCE_DIR}/../lib/cPMML cPMML)

# Create an external library target
#add_library(cPMML STATIC IMPORTED)
#
## Set the library location (will be set after build)
#set_target_properties(cPMML PROPERTIES
#        IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/src/cPMML/libcPMML.a
#)
## Ensure the library is built before linking
#add_dependencies(cPMML protal_lib)

#find_package(cPMML REQUIRED)


#set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")  # Set to OFF for static build


#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#find_library(CPMML_LIBRARY
#        NAMES libcPMML_lib.a
#        PATHS
#        ${CMAKE_BINARY_DIR}
#        NO_DEFAULT_PATH
#)

#set(CPMML_LIBRARY "${CMAKE_BINARY_DIR}/src/cPMML/libcPMML.a")



message("ProjectSource: ${PROJECT_SOURCE_DIR}")
message("CPMML_LIBRARY: ${CPMML_LIBRARY}")
message("CMAKE PREFIX PATH: ${CMAKE_PREFIX_PATH}")
message("Find cPMML: ${PROJECT_SOURCE_DIR}/../")
message("Build dir: ${CMAKE_BINARY_DIR}/src/cPMML/")
message("CPMML: ${cPMML}")


# Add this before target_link_libraries
find_library(RT_LIBRARY rt)

message("RT LIB: ${RT_LIBRARY}")

# Necessary to link to stdc++fs (filesystem experimental support)
target_link_libraries(protal_lib
        LINK_PUBLIC
        ${RT_LIBRARY}
        cPMML
        stdc++fs
        wfa_lib
        gzstream_lib
        ZLIB::ZLIB
        -static # WITHOUT THIS IT WILL NOT COMPILE INTO A STATIC LIBRARY
)

message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


message("Ranlib: ${CMAKE_RANLIB}")
message("Ranlib: ${CMAKE_AR}")
message("GCC: $ENV{GCC}")
