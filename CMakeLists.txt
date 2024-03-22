# minimum required CMake version
cmake_minimum_required(VERSION 3.18.0)

# if on a Mac, use gnu compilers instead of clang
if(APPLE)
    set(CMAKE_C_COMPILER gcc-13)
    set(CMAKE_CXX_COMPILER g++-13)
endif()

# set the project name and get version from version.txt
file(READ "grss/version.txt" ver)
message(STATUS "GRSS version: ${ver}")
project(grss VERSION ${ver})

# specify the C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-std=c++11 -O3 -fPIC -fopenmp) # operational flags
# add_compile_options(-std=c++11 -g3 -fPIC -fopenmp -Werror -Wall -Wextra -pedantic) # debugging flags

# Set header file directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/cspice/include)

# Set source code directories
add_subdirectory(src)

# create pybind11 module
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)
pybind11_add_module(libgrss src/${PROJECT_NAME}.cpp)
target_link_libraries(libgrss PRIVATE ${PROJECT_NAME})
target_compile_definitions(libgrss PRIVATE VERSION_INFO=${${PROJECT_NAME}_VERSION})

# compile tests
add_subdirectory(tests/cpp/prop)
