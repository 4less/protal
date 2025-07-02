cmake_minimum_required(VERSION 3.22)
project(protal VERSION 0.1.0)

option(BUILD_BENCHMARKS "Build benchmarks" OFF)
if(BUILD_BENCHMARKS)
  add_executable(model_benchmark.exe ...)
  # ... rest of benchmark setup ...
endif()

message("CONDA_PREFIX: $ENV{CONDA_PREFIX}")
message("PATH: $ENV{PATH}")
message("HOST: $ENV{HOST}")
message("PREFIX: $ENV{PREFIX}")
message("PROJECT_SOURCE: ${PROJECT_SOURCE_DIR}")
message("BUILD_BENCHMARKS: ${BUILD_BENCHMARKS}")


# god-awful hack because it seems to not run correct tests to determine this:
set(__CHAR_UNSIGNED___EXITCODE 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_FLAGS}")

#message("CPLUS_INCLUDE_PATH")
#message("$ENV{CPLUS_INCLUDE_PATH}")


configure_file(protal_config.h.in ${PROJECT_SOURCE_DIR}/src/protal_config.h)

message("CONDA_PREFIX: $ENV{CONDA_PREFIX}")
message("PATH: $ENV{PATH}")
message("HOST: $ENV{HOST}")
message("PREFIX: $ENV{PREFIX}")
message("ProjectSource: ${PROJECT_SOURCE_DIR}/src/protal_config.h")

set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_EXE_LINKER_FLAGS "-lz")
# set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a")# -mavx -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++2a -DNDEBUG -pthread -ggdb -g -fopenmp -Werror=unknown-pragmas") #-march=native
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native


#set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -mavx2 -fopenmp")
#set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -mavx2 -fopenmp")

include_directories(src/ lib/)

add_executable(protal
        src/main.cpp
        protal_config.h.in)

# Tests

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(test_CigarIterator tests/test_CigarIterator.cpp)
target_link_libraries(test_CigarIterator GTest::gtest_main)
target_include_directories(test_CigarIterator PRIVATE src)
gtest_discover_tests(test_CigarIterator)

add_executable(test_Coverage tests/test_Coverage.cpp)
target_link_libraries(test_Coverage GTest::gtest_main protal_lib)
target_include_directories(test_Coverage PRIVATE src)
gtest_discover_tests(test_Coverage)
# Tests end

set_target_properties(protal PROPERTIES COMPILE_FLAGS "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas")

add_executable(protal_avx2
        src/main.cpp
        protal_config.h.in)
#target_compile_options(protal_avx2 PRIVATE -std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx)
set_target_properties(protal_avx2 PROPERTIES COMPILE_FLAGS "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx")




# Link to cmake file in src
target_link_libraries(protal PRIVATE protal_lib)
target_link_libraries(protal_avx2 PRIVATE protal_lib)

add_subdirectory(src/)
#add_subdirectory(tests/)
