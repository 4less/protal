#message("GOOFER?")
#message("GOOFBALL: $ENV{GOOFBALL}")
#message("CONDA_PREFIX: $ENV{CONDA_PREFIX}")
#message("PATH: $ENV{PATH}")
#message("HOST: $ENV{HOST}")
#message("PREFIX: $ENV{PREFIX}")

# where is the target environment
#set(CMAKE_FIND_ROOT_PATH $ENV{PREFIX} $ENV{CONDA_PREFIX}/$ENV{HOST}/sysroot)
#message("CMAKE FIND ROOT PATH: ${CMAKE_FIND_ROOT_PATH}")
#
# search for programs in the build host directories
#set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
#set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


project(gzstream_lib)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXE_LINKER_FLAGS "-lgzstream -lz")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native

set(ZLIB_USE_STATIC_LIBS "ON")
find_package(ZLIB REQUIRED)

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
        PUBLIC
        ZLIB::ZLIB)