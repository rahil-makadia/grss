# minimum required CMake version
cmake_minimum_required(VERSION 3.19.0)

# set C and C++ compilers for Apple M1 - openmp issues
# if (APPLE)
#     set(CMAKE_C_COMPILER "gcc-12")
#     set(CMAKE_CXX_COMPILER "g++-12")
# endif()

# set the project name and version
project(grss VERSION 0.1.0)

# specify the C++ standard and compiler flags
# TODO: add gprof profiling, iee standards for floats
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-std=c++11 -DLONGDOUBLE -O3 -fPIC -Werror -Wall -Wextra -pedantic -Wno-unused-but-set-variable)

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
pybind11_add_module(cpp${PROJECT_NAME} src/${PROJECT_NAME}.cpp)
target_link_libraries(cpp${PROJECT_NAME} PRIVATE ${PROJECT_NAME})
target_compile_definitions(cpp${PROJECT_NAME} PRIVATE VERSION_INFO=${${PROJECT_NAME}_VERSION})

# copy pybind11 module to grss directory
add_custom_command(TARGET cpp${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:cpp${PROJECT_NAME}> ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/$<TARGET_FILE_NAME:cpp${PROJECT_NAME}>
    # COMMENT "Copying pybind11 module to grss directory"
)
