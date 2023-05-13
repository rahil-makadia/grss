add_library(${PROJECT_NAME}
force.cpp
gr15.cpp
interpolate.cpp
simulation.cpp
spk.cpp
utilities.cpp
)

# link cspice
if (WIN32)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cspice/lib/cspice.lib)
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cspice/lib/cspice.a)
endif()

# link openmp
# find_package(OpenMP REQUIRED)
# if(OPENMP_FOUND)
#     target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
# endif()
