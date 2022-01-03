set(SRCS
    pyclasses/bspline1.cpp
    pyclasses/bspline2.cpp
    pyclasses/bspline3.cpp
    pyclasses/bspline4.cpp
    pyclasses/bspline5.cpp
    pyclasses/bspline6.cpp
    pyclasses/bspline7.cpp
    pyclasses/bspline8.cpp
    pyclasses/bspline9.cpp
    pyclasses/bspline10.cpp
    pyclasses/nurbs1.cpp
    pyclasses/nurbs2.cpp
    pyclasses/nurbs3.cpp
    pyclasses/nurbs4.cpp
    pyclasses/nurbs5.cpp
    pyclasses/nurbs6.cpp
    pyclasses/nurbs7.cpp
    pyclasses/nurbs8.cpp
    pyclasses/nurbs9.cpp
    pyclasses/nurbs10.cpp
    pyclasses/reader.cpp
    splinepy.cpp
)

#set(SRCS splinelib.cpp)

pybind11_add_module(_splinepy ${SRCS})
target_include_directories(_splinepy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

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
endif()

# Link Fitting
add_subdirectory(fitting)
target_link_libraries(_splinepy PRIVATE fitting)
