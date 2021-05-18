cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU C CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS regex graph)

add_subdirectory(lib/mmio)
add_subdirectory(lib/zf_log)
include_directories(lib/snap-core)
include_directories(lib/glib-core)
include_directories(lib/doctest)
include_directories(include)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(dataset)

#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")
#set(CMAKE_C_FLAGS_RELEASE "-O2 -DBENCHMARK")

# overwrote on purpose
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g -G -gencode arch=compute_60,code=sm_60 -Wall -Wpedantic -Wextra")

#message("-- CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
#message("-- CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")

message("-- CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("-- CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    add_compile_definitions(try -DNDEBUG -DRELEASE -DNO_CUDA_CHECK_ERROR)
endif ()
