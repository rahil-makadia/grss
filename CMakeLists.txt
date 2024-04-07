# minimum required CMake version
cmake_minimum_required(VERSION 3.18.0)

# if on a Mac, use gnu compilers instead of clang for ease of use with OpenMP
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
add_compile_options(-std=c++11 -O3 -fPIC -fopenmp -DGRSS_VERSION="${${PROJECT_NAME}_VERSION}") # operational flags
# add_compile_options(-std=c++11 -g3 -fPIC -fopenmp -DGRSS_VERSION="${${PROJECT_NAME}_VERSION}" -Werror -Wall -Wextra -pedantic) # debugging flags

# Set header file directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Set source code directory
add_subdirectory(src)

# create pybind11 module
find_package(pybind11 REQUIRED)
pybind11_add_module(libgrss src/${PROJECT_NAME}.cpp)
target_link_libraries(libgrss PRIVATE ${PROJECT_NAME})

# compile tests
add_subdirectory(tests/cpp/prop)

# build the documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/doxygen/doxygen.config)
    set(DOXYGEN_OUT ${CMAKE_SOURCE_DIR}/docs/doxygen/doxyfile.out)

    # configure the file and fill out CMake variables
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # note: the option ALL allows building the docs with the binaries
    add_custom_target(
        doxygen_docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs/doxygen
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
else (DOXYGEN_FOUND)
    message(WARNING "Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
