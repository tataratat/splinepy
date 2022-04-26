# spline sources
if(MINIMAL)
    set(SRCS python/minimal.cpp)
    set(DEFS _MINIMAL_)

else()
    set(
        SRCS
            python/bezier1.cpp
            python/bezier2.cpp
            python/bezier3.cpp
            python/bezier4.cpp
            python/bezier5.cpp
            python/bezier6.cpp
            python/bezier7.cpp
            python/bezier8.cpp
            python/bezier9.cpp
            python/bezier10.cpp
            python/bspline1.cpp
            python/bspline2.cpp
            python/bspline3.cpp
            python/bspline4.cpp
            python/bspline5.cpp
            python/bspline6.cpp
            python/bspline7.cpp
            python/bspline8.cpp
            python/bspline9.cpp
            python/bspline10.cpp
            python/nurbs1.cpp
            python/nurbs2.cpp
            python/nurbs3.cpp
            python/nurbs4.cpp
            python/nurbs5.cpp
            python/nurbs6.cpp
            python/nurbs7.cpp
            python/nurbs8.cpp
            python/nurbs9.cpp
            python/nurbs10.cpp
    )

endif(MINIMAL)

# reader and python module
list(APPEND SRCS python/reader.cpp splinepy.cpp)


pybind11_add_module(_splinepy ${SRCS})
target_include_directories(_splinepy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(_splinepy PRIVATE ${DEFS})

target_compile_options(_splinepy PRIVATE ${SPLINELIBPY_FLAGS})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(_splinepy PRIVATE ${GNU_FLAGS_R})
endif()

# Link SplineLib
if(COMPILE_SPLINELIB)
    target_link_libraries(_splinepy PRIVATE splines input_output)
else()
    find_package(SplineLib REQUIRED)
    target_link_libraries(_splinepy 
        PRIVATE
            SplineLib::splines
            SplineLib::input_output
    )
    add_library(splines ALIAS SplineLib::splines)
    add_library(input_output ALIAS SplineLib::input_output)
endif()

# Link bezierManipulation
target_link_libraries(_splinepy PRIVATE bezierManipulation)

# Link Fitting
add_subdirectory(fitting)
target_link_libraries(_splinepy PRIVATE fitting)
