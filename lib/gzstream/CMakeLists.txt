project(gzstream_lib)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXE_LINKER_FLAGS "-lgzstream -lz")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native

file (GLOB gzstream_sources
        "*.c"
        "*.h"
        "*.C")

message("GZSTREAM-Sources")
message("${gzstream_sources}")

# This creates the target and must be specified before using target_include etc
add_library(gzstream_lib
        ${gzstream_sources})

set_target_properties(gzstream_lib PROPERTIES LINKER_LANGUAGE CXX)

target_link_directories(gzstream_lib
        PUBLIC
        ${PROJECT_SOURCE_DIR})

target_include_directories(gzstream_lib
        PUBLIC
        ${PROJECT_SOURCE_DIR}/)

target_link_libraries(gzstream_lib
        LINK_PUBLIC
        ZLIB::ZLIB)