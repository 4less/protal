cmake_minimum_required(VERSION 3.16)
project(wfa2lib)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXE_LINKER_FLAGS "-lz")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -march=native -mavx2")


file (GLOB wfa_sources
        "alignment/*.c"
        "alignment/*.h"
        "utils/*.h"
        "utils/*.c"
        bindings/cpp/WFAligner.hpp
        bindings/cpp/WFAligner.cpp
        scripts/fasta2seq.sh
        scripts/seq2fasta.sh
        "system/*.c"
        "system/*.h"
        "wavefront/*.h"
        "wavefront/*.c"
)

message("WFA-Sources")
message("${wfa_sources}")

# This creates the target and must be specified before using target_include etc
add_library(wfa2lib
        ${wfa_sources}
        #        alignment/
        #        wavefront/
        #        system/
        #        scripts/
        #        bindings/
        #        utils/
)


set_target_properties(wfa2lib PROPERTIES LINKER_LANGUAGE CXX)


###########################################################################################
# Now target varkit_lib exists
###########################################################################################

# link_directories makes files in these directories includable without explicitly stating the dir
# include "test/kmer_iterator_test.h" -> include "kmer_iterator_test.h"

message(${PROJECT_SOURCE_DIR}/)

target_link_directories(wfa2lib
        PUBLIC
        ${PROJECT_SOURCE_DIR}/
        wavefront/
        alignment/
        system/
        scripts/
        bindings/
        bindings/cpp/
        utils/)

target_include_directories(wfa2lib
        PUBLIC
        ${PROJECT_SOURCE_DIR}/
        wavefront/
        alignment/
        system/
        scripts/
        bindings/
        bindings/cpp/
        utils/)
