# disable tests for splinelib and bezman
if(SPLINEPY_COMPILE_BSPLINELIB)
  add_subdirectory(BSplineLib)
else()
  # check here if bsplinelib is available from cmake 3.24 we could add `GLOBAL`
  # flag
  find_package(BSplineLib REQUIRED)
endif()

set(BEZMAN_GOOGLETEST OFF)
set(BEZMAN_USE_COMPILE_OPTIONS OFF)
add_subdirectory(bezman)

add_subdirectory(pybind11)
add_subdirectory(napf)
