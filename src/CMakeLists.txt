add_library(${PROJECT_NAME}
approach.cpp
elements.cpp
force.cpp
ias15.cpp
interpolate.cpp
observe.cpp
parallel.cpp
pck.cpp
simulation.cpp
spk.cpp
stm.cpp
timeconvert.cpp
utilities.cpp
)

# link openmp
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
