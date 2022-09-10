cmake_minimum_required(VERSION 3.22)
project(protal VERSION 0.0.1)


configure_file(protal_config.h.in ${PROJECT_SOURCE_DIR}/src/protal_config.h)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXE_LINKER_FLAGS "-lz")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native

include_directories(src/ lib/)

add_executable(protal
        src/main.cpp
        protal_config.h.in)

# Link to cmake file in src
target_link_libraries(protal protal_lib)

add_subdirectory(src/)
#add_subdirectory(tests/)