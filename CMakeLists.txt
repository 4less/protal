cmake_minimum_required(VERSION 3.22)
project(protal VERSION 0.6.0)

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
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -Werror=unknown-pragmas") #-march=native
set(OPENMP_FLAG "")


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
	if (CMAKE_DISABLE_FIND_PACKAGE_OpenMP)
		message(STATUS "OpenMP disabled by configuration")
	else()
		find_package(OpenMP REQUIRED)
		set (OPENMP_FLAG "-fopenmp")
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ${OPENMP_FLAG}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${OPENMP_FLAG}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

if (CMAKE_DISABLE_FIND_PACKAGE_OpenMP)
	# Avoid failing builds when OpenMP pragmas are present but support is disabled.
	add_compile_options(-Wno-error=unknown-pragmas)
	string(REPLACE "-Werror=unknown-pragmas" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
	string(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(PROTAL_COMMON_COMPILE_OPTIONS
        -std=c++2a
        -DNDEBUG
        -pthread
        -O3
        -Werror=unknown-pragmas)
if (OPENMP_FLAG)
	list(APPEND PROTAL_COMMON_COMPILE_OPTIONS ${OPENMP_FLAG})
endif()

add_library(isa_baseline INTERFACE)
target_compile_options(isa_baseline INTERFACE -march=x86-64 -mtune=generic)

add_library(isa_avx2 INTERFACE)
target_compile_options(isa_avx2 INTERFACE -march=x86-64-v3 -mtune=generic)

include_directories(src/ lib/)

add_executable(protal
        src/main.cpp
        protal_config.h.in)
target_compile_options(protal PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(protal_static
        src/main.cpp
        protal_config.h.in)
target_compile_options(protal_static PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_options(protal_static PRIVATE -static -static-libgcc -static-libstdc++)

add_executable(protal_avx2
        src/main.cpp
        protal_config.h.in)
target_compile_options(protal_avx2 PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(simulate_metagenomes
        src/simulate_metagenomes_main.cpp)
target_compile_options(simulate_metagenomes PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(simulate_metagenomes_static
        src/simulate_metagenomes_main.cpp)
target_compile_options(simulate_metagenomes_static PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_options(simulate_metagenomes_static PRIVATE -static -static-libgcc -static-libstdc++)


# Link to cmake file in src
target_link_libraries(protal PRIVATE protal_lib dl z isa_baseline)
target_link_libraries(protal_static PRIVATE protal_lib_static dl z isa_baseline)
target_link_libraries(protal_avx2 PRIVATE protal_lib dl z isa_avx2)
target_link_libraries(simulate_metagenomes PRIVATE protal_lib dl z isa_baseline)
target_link_libraries(simulate_metagenomes_static PRIVATE protal_lib_static dl z isa_baseline)

add_subdirectory(src/)
#add_subdirectory(tests/)
