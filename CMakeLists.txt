# minimum required CMake version
cmake_minimum_required(VERSION 3.19.0)

# set C and C++ compilers for Apple M1 - openmp issues
if (APPLE)
    set(CMAKE_C_COMPILER "gcc-12")
    set(CMAKE_CXX_COMPILER "g++-12")
endif()

# set the project name and version
project(grss VERSION 0.1.0)

# specify the C++ standard and compiler flags
# TODO: add gprof profiling, iee standards for floats
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-std=c++11 -DLONGDOUBLE -O3 -fPIC -Wall -Wextra -Werror -pedantic) # -Weffc++

# Set output directories, from https://stackoverflow.com/a/6595001
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Set header file directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/cspice/include)
# include_directories(${CMAKE_SOURCE_DIR}/extern/pybind11/include)

# Set source code directories
add_subdirectory(extern/pybind11)
add_subdirectory(src)
add_subdirectory(examples)

# Configure CPack for cross-platform packaging
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# create pybind11 module
pybind11_add_module(${PROJECT_NAME}_py src/${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME}_py PRIVATE ${PROJECT_NAME})
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a define (VERSION_INFO) here.
target_compile_definitions(${PROJECT_NAME}_py PRIVATE VERSION_INFO=${${PROJECT_NAME}_VERSION})
