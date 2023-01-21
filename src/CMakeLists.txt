project(protal_lib)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXE_LINKER_FLAGS " -lwfa -lz")
set(CMAKE_CXX_FLAGS_DEBUG " -pthread -ggdb -g -std=c++2a -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++20 -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native

message(${CMAKE_CXX_FLAGS_DEBUG})

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
        Profiling/ProfilerOptions.h
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
        )


set_target_properties(protal_lib PROPERTIES LINKER_LANGUAGE CXX)

function(get_all_targets var)
    set(targets)
    get_all_targets_recursive(targets ${CMAKE_CURRENT_SOURCE_DIR})
    set(${var} ${targets} PARENT_SCOPE)
endfunction()

macro(get_all_targets_recursive targets dir)
    get_property(subdirectories DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    foreach(subdir ${subdirectories})
        get_all_targets_recursive(${targets} ${subdir})
    endforeach()

    get_property(current_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    list(APPEND ${targets} ${current_targets})
endmacro()

get_all_targets(all_targets)
message("All targets: ${all_targets}")


###########################################################################################
# Now target varkit_lib exists
###########################################################################################

# link_directories makes files in these directories includable without explicitly stating the dir
# include "test/kmer_iterator_test.h" -> include "kmer_iterator_test.h"

message(${PROJECT_SOURCE_DIR})

target_link_directories(protal_lib
        PUBLIC
        ../lib/
        IO/
        Hash/
        Core/
        VarkitInterface/
        Alignment/
        Utils/
        Profiling/
        SequenceUtils/
        Taxonomy/
        SNPs/)

target_include_directories(protal_lib
        PUBLIC
        ../lib/
        IO/
        Hash/
        Core/
        VarkitInterface/
        Alignment/
        Utils/
        Profiling/
        SequenceUtils/
        Taxonomy/
        SNPs/)


# Adding this to the path makes find_package find the required BHM library
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
#link_directories(/usr/users/QIB_fr017/fritsche/library_cpp/bhm/lib)


#message(STATUS "BHM include dir: ${BHM_INCLUDE_DIR}")
#message(STATUS "BHM library: ${BHM_LIBRARY}")
#message(STATUS "CSV library: ${CWCSV_INCLUDE_DIR}")

# Directly include folders to make them "includable" with '#include "robin_map.h"'
# require BitHashMap
# find_package(BHM REQUIRED)

find_package(ZLIB REQUIRED)
#find_library(HTSLIB htslib /usr/users/QIB_fr017/fritsche/miniconda3/include/)

find_library(HTSLIB_LIB hts)
target_link_libraries(protal_lib PUBLIC ${HTSLIB_LIB})

# include CWCSV
#find_package(CWCSV REQUIRED)

# Directly include folders to make them "includable" with '#include "robin_map.h"'
# generic include_directories should link with all targets in "dir"
# however this did not work, so I had to change to target_include_directories

# Both needed (for whatever reason)
#target_include_directories(varkit_lib PUBLIC ${BHM_INCLUDE_DIR})
#target_link_libraries(varkit_lib LINK_PUBLIC ${BHM_LIBRARY})

target_include_directories(protal_lib PUBLIC /usr/users/QIB_fr017/fritsche/library_cpp/robin)
target_include_directories(protal_lib PUBLIC /usr/users/QIB_fr017/fritsche/library_cpp/tsl)


# Necessary to link to stdc++fs (filesystem experimental support)
target_link_libraries(protal_lib
        LINK_PUBLIC
        ZLIB::ZLIB
        stdc++fs
        wfa_lib
        gzstream_lib
#        hts_lib
        #fmt_lib
        )

add_subdirectory(/usr/users/QIB_fr017/fritsche/CLionProjects/protal/lib/wfa2-lib/ /usr/users/QIB_fr017/fritsche/CLionProjects/protal/src/wfa)
#add_subdirectory(/usr/users/QIB_fr017/fritsche/CLionProjects/protal/lib/fmt/ /usr/users/QIB_fr017/fritsche/CLionProjects/protal/src/fmt)
add_subdirectory(/usr/users/QIB_fr017/fritsche/CLionProjects/protal/lib/gzstream/ /usr/users/QIB_fr017/fritsche/CLionProjects/protal/src/gzstream)
#add_subdirectory(/usr/users/QIB_fr017/fritsche/CLionProjects/protal/lib/htslib-1.16/ /usr/users/QIB_fr017/fritsche/CLionProjects/protal/src/htslib)


#target_include_directories(protal_lib PUBLIC ${CWCSV_INCLUDE_DIR})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CWCSV_INCLUDE_DIR})