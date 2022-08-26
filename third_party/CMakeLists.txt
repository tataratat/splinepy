# disable tests for splinelib and bezman
set(GOOGLETEST OFF)

if(SPLINEPY_COMPILE_SPLINELIB)
    add_subdirectory(SplineLib)
endif()
add_subdirectory(bezman)
add_subdirectory(pybind11)
add_subdirectory(napf)
