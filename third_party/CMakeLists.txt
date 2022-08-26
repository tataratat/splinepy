# disable tests for splinelib and bezman
set(GOOGLETEST OFF)

if(SPLINEPY_COMPILE_SPLINELIB)
    add_subdirectory(SplineLib)
endif()
#add_subdirectory(bezman)
# temporary workaround to avoid unwanted flag propagation
add_library(bezman INTERFACE)
target_include_directories(bezman INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(pybind11)
add_subdirectory(napf)
