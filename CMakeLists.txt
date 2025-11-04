cmake_minimum_required(VERSION 3.22)
project(protal VERSION 0.2.1)

message("CONDA_PREFIX: $ENV{CONDA_PREFIX}")
message("PATH: $ENV{PATH}")
message("HOST: $ENV{HOST}")
message("PREFIX: $ENV{PREFIX}")
message("PROJECT_SOURCE: ${PROJECT_SOURCE_DIR}")


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
#set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a")# -mavx -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native


#set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -mavx2 -fopenmp")
#set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -mavx2 -fopenmp")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set (OpenMP_C "${CMAKE_C_COMPILER}")
	set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
	set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
	set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
	set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
	set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
	set (OpenMP_libomp_LIBRARY "omp")
	set (OpenMP_libgomp_LIBRARY "gomp")
	set (OpenMP_libiomp5_LIBRARY "iomp5")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lomp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	find_package(OpenMP REQUIRED)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(src/ lib/)

add_executable(protal
        src/main.cpp
        protal_config.h.in)
set_target_properties(protal PROPERTIES COMPILE_FLAGS "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas")

add_executable(protal_avx2
        src/main.cpp
        protal_config.h.in)
#target_compile_options(protal_avx2 PRIVATE -std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx)
set_target_properties(protal_avx2 PROPERTIES COMPILE_FLAGS "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas -mavx2 -mavx")



# Link to cmake file in src
target_link_libraries(protal PRIVATE protal_lib dl z)
target_link_libraries(protal_avx2 PRIVATE protal_lib dl z)

add_subdirectory(src/)
#add_subdirectory(tests/)

