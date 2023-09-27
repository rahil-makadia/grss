# minimum required CMake version
cmake_minimum_required(VERSION 3.18.0)

# set the project name and get version from version.txt
file(READ "grss/version.txt" ver)
message(STATUS "GRSS version: ${ver}")
project(grss VERSION ${ver})

# specify the C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-std=c++11 -O3 -fPIC) # operational flags
# add_compile_options(-std=c++11 -DLONGDOUBLE -g3 -fPIC -Werror -Wall -Wextra -pedantic -Wno-unused-but-set-variable -Wno-unused-result) # debugging flags

# Set header file directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/cspice/include)

# Set source code directories
add_subdirectory(extern/pybind11)
add_subdirectory(src)

# Configure CPack for cross-platform packaging
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# create pybind11 module
pybind11_add_module(prop_simulation src/${PROJECT_NAME}.cpp)
target_link_libraries(prop_simulation PRIVATE ${PROJECT_NAME})
target_compile_definitions(prop_simulation PRIVATE VERSION_INFO=${${PROJECT_NAME}_VERSION})

# compile tests
add_subdirectory(tests/cpp/prop)
