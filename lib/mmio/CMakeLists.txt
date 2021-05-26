cmake_minimum_required(VERSION 3.2)

set(HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(HEADERS mmio.h)
set(SOURCES mmio.c)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(mmio ${HEADERS} ${SOURCES})
target_include_directories(mmio PUBLIC $<BUILD_INTERFACE:${HEADERS_DIR}>)
