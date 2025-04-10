# minimum required CMake version
cmake_minimum_required(VERSION 3.18.0)
# if cmake version is greater than 3.27, set the policy to OLD
if(${CMAKE_VERSION} VERSION_GREATER "3.27")
    cmake_policy(SET CMP0148 OLD)
endif()

# if on a Mac, use any GNU compiler from Homebrew
if(APPLE)
    # find the brew GNU compiler version
    execute_process(COMMAND bash -c "ls $HOMEBREW_PREFIX/bin/g++*" OUTPUT_VARIABLE BREW_GXX OUTPUT_STRIP_TRAILING_WHITESPACE)
    # last two characters are the version number
    string(REGEX REPLACE ".*g\\+\\+-(.*)$" "\\1" BREW_GXX_VERSION ${BREW_GXX})
    message(STATUS "Found GNU g++ compiler version: ${BREW_GXX_VERSION}")
    # use any GNU g++ compiler in BREW_PREFIX/bin without specifying the version
    set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/bin/gcc-${BREW_GXX_VERSION}")
    set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/bin/g++-${BREW_GXX_VERSION}")
endif()

# set the project name and get version from version.txt
file(READ "grss/version.txt" ver)
message(STATUS "GRSS version: ${ver}")
project(grss VERSION ${ver})

# specify the C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-std=c++11 -O3 -fPIC -fopenmp -DGRSS_VERSION="${${PROJECT_NAME}_VERSION}") # operational flags
# add_compile_options(-std=c++11 -g2 -fPIC -fopenmp -DGRSS_VERSION="${${PROJECT_NAME}_VERSION}" -Werror -Wall -Wextra -pedantic) # debugging flags

# Set static library output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# Set header file directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Set source code directory
add_subdirectory(src)

# create pybind11 module
find_package(pybind11 REQUIRED)
pybind11_add_module(libgrss src/${PROJECT_NAME}.cpp)
target_link_libraries(libgrss PRIVATE ${PROJECT_NAME})

# copy the libgrss* files to the python package directory
file(GLOB LIBGRSS_FILES ${CMAKE_SOURCE_DIR}/build/libgrss.*)
file(COPY ${LIBGRSS_FILES} DESTINATION ${CMAKE_SOURCE_DIR}/grss)

# compile tests
add_subdirectory(tests/cpp/prop)

# build the documentation if Doxygen is installed and BUILD_DOCS is defined
find_package(Doxygen)
if (DOXYGEN_FOUND AND BUILD_DOCS)
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
else (DOXYGEN_FOUND AND BUILD_DOCS)
    message(STATUS "Doxygen needs to be installed and BUILD_DOCS needs to be ON to build the documentation")
endif (DOXYGEN_FOUND AND BUILD_DOCS)
