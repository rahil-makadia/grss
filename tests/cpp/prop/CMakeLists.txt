# create executables for example files and link against grss library
file(GLOB files "${CMAKE_SOURCE_DIR}/tests/cpp/prop/*.cpp")
foreach(file ${files})
    get_filename_component(name ${file} NAME_WE)
    add_executable(${name} ${file})
    set_target_properties(${name} PROPERTIES SUFFIX ".out")
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/cpp/prop)
    target_link_libraries(${name} ${PROJECT_NAME})
endforeach()
