add_library(${PROJECT_NAME}
approach.cpp
elements.cpp
force.cpp
gr15.cpp
interpolate.cpp
parallel.cpp
pck.cpp
simulation.cpp
spk.cpp
stm.cpp
timeconvert.cpp
utilities.cpp
)

# link cspice
if (WIN32)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cspice/lib/cspice.lib)
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/extern/cspice/lib/cspice.a)
endif()

# link openmp
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
