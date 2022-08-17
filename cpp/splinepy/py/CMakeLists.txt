if(MINIMAL)
    set(SRCS init/minimal.cpp)
    set(DEFS _MINIMAL_)

else()
    set(
        SRCS
            init/py_bezier1.cpp
            init/py_bezier2.cpp
            init/py_bezier3.cpp
            init/py_bezier4.cpp
            init/py_bezier5.cpp
            init/py_bezier6.cpp
            init/py_bezier7.cpp
            init/py_bezier8.cpp
            init/py_bezier9.cpp
            init/py_bezier10.cpp
            init/py_rational_bezier1.cpp
            init/py_rational_bezier2.cpp
            init/py_rational_bezier3.cpp
            init/py_rational_bezier4.cpp
            init/py_rational_bezier5.cpp
            init/py_rational_bezier6.cpp
            init/py_rational_bezier7.cpp
            init/py_rational_bezier8.cpp
            init/py_rational_bezier9.cpp
            init/py_rational_bezier10.cpp
            init/py_bspline1.cpp
            init/py_bspline2.cpp
            init/py_bspline3.cpp
            init/py_bspline4.cpp
            init/py_bspline5.cpp
            init/py_bspline6.cpp
            init/py_bspline7.cpp
            init/py_bspline8.cpp
            init/py_bspline9.cpp
            init/py_bspline10.cpp
            init/py_nurbs1.cpp
            init/py_nurbs2.cpp
            init/py_nurbs3.cpp
            init/py_nurbs4.cpp
            init/py_nurbs5.cpp
            init/py_nurbs6.cpp
            init/py_nurbs7.cpp
            init/py_nurbs8.cpp
            init/py_nurbs9.cpp
            init/py_nurbs10.cpp
    )

endif(MINIMAL)

# reader and pyhton module
list(APPEND SRCS init/reader.cpp init/exporter.cpp splinepy.cpp)

pybind11_add_module(_splinepy ${SRCS})
target_include_directories(_splinepy PRIVATE "${PROJECT_SOURCE_DIR}/cpp")
target_compile_definitions(_splinepy PRIVATE ${DEFS})

target_compile_options(_splinepy PRIVATE
    ${SPLINEPY_FLAGS}
    ${SPLINEPY_OPTIMIZATION_FLAGS})

# link splinelib
if(SPLINEPY_COMPILE_SPLINELIB)
    target_link_libraries(_splinepy PRIVATE splines input_output)
else()
    find_package(SplineLib REQUIRED)
    target_link_libraries(_splinepy
        PRIVATE
            SplineLib::splines
            SplineLib::input_output
    )
endif()

# link bezman
target_link_libraries(_splinepy PRIVATE bezman)

# link fitting
target_link_libraries(_splinepy PRIVATE fitting)

# link napf
target_link_libraries(_splinepy PRIVATE napf)
