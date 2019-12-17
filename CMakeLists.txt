cmake_minimum_required(VERSION 2.8.8)

set(PROJECT_NAME_STR smart)
PROJECT(${PROJECT_NAME_STR} C CXX)
#TODO CHANGE THIS OPTION
set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_BUILD_TYPE Release)

find_package(Threads REQUIRED)

execute_process(COMMAND ${PROJECT_SOURCE_DIR}/pre-install.sh)

# SDSL directory
set(SDSL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/ext/sdsl-lite")

# flags
set(CMAKE_CXX_FLAGS               "-Wall -Wextra -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# include the sdsl files and googletest library
add_subdirectory("${PROJECT_SOURCE_DIR}/ext/libdivsufsort")
add_subdirectory("${PROJECT_SOURCE_DIR}/ext/sdsl-lite")

# add libdivsufsort includes path
include_directories("${CMAKE_BINARY_DIR}/ext/libdivsufsort/include")
include_directories("${CMAKE_BINARY_DIR}/ext/sdsl-lite/include/sdsl")

# Compile project
ADD_EXECUTABLE(${PROJECT_NAME_STR}
               src/sa.cpp
               src/ReadsDB.cpp src/AppConfig.cpp
               src/ExactLCPk.cpp
               src/compute_klcp.cpp
               src/main.cpp src/matrices.cpp)

#Link SDSL library
target_link_libraries(${PROJECT_NAME_STR} divsufsort divsufsort64)
target_link_libraries(${PROJECT_NAME_STR} sdsl)
