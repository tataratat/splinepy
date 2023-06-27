# disable tests for splinelib and bezman
set(GOOGLETEST OFF)

if(SPLINEPY_COMPILE_SPLINELIB)
    add_subdirectory(SplineLib)
endif()

set(BEZMAN_GOOGLETEST OFF)
set(BEZMAN_USE_COMPILE_OPTIONS OFF)
add_subdirectory(bezman)

add_subdirectory(pybind11)
add_subdirectory(napf)
