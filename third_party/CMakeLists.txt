# disable tests for splinelib and bezman
if(SPLINEPY_COMPILE_BSPLINELIB)
    add_subdirectory(BSplineLib)
endif()

set(BEZMAN_GOOGLETEST OFF)
set(BEZMAN_USE_COMPILE_OPTIONS OFF)
add_subdirectory(bezman)

add_subdirectory(pybind11)
add_subdirectory(napf)
