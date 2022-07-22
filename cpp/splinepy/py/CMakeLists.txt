if(MINIMAL)
    set(SRCS init/minimal.cpp)
    set(DEFS _MINIMAL_)

else()
    set(
        SRCS
            init/bezier1.cpp
            init/bezier2.cpp
            init/bezier3.cpp
            init/bezier4.cpp
            init/bezier5.cpp
            init/bezier6.cpp
            init/bezier7.cpp
            init/bezier8.cpp
            init/bezier9.cpp
            init/bezier10.cpp
            init/bspline1.cpp
            init/bspline2.cpp
            init/bspline3.cpp
            init/bspline4.cpp
            init/bspline5.cpp
            init/bspline6.cpp
            init/bspline7.cpp
            init/bspline8.cpp
            init/bspline9.cpp
            init/bspline10.cpp
            init/nurbs1.cpp
            init/nurbs2.cpp
            init/nurbs3.cpp
            init/nurbs4.cpp
            init/nurbs5.cpp
            init/nurbs6.cpp
            init/nurbs7.cpp
            init/nurbs8.cpp
            init/nurbs9.cpp
            init/nurbs10.cpp
    )

endif(MINIMAL)

# reader and pyhton module
list(APPEND SRCS init/reader.cpp splinepy.cpp)

pybind11_add_module(_splinepy ${SRCS})
target_include_directories(_splinepy PRIVATE "${PROJECT_SOURCE_DIR}/cpp")
target_compile_definitions(_splinepy PRIVATE ${DEFS})

target_compile_options(_splinepy PRIVATE
    ${SPLINEPY_FLAGS}
    ${SPLINEPY_OPTIMIZATION_FLAGS})

# link splinelib
if(COMPILE_SPLINELIB)
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
target_link_libraries(_splinepy PRIVATE bezierManipulation)

# link fitting
target_link_libraries(_splinepy PRIVATE fitting)

# link napf
target_link_libraries(_splinepy PRIVATE napf)
