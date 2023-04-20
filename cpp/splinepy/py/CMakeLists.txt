# reader, export helper, and pyhton module
set(PYSPLINEPY_SRCS
        init/reader.cpp
        init/exporter.cpp
        init/coordinate_reference.cpp
        init/core_spline.cpp
        init/fitting.cpp
        init/spline_extensions.cpp
        splinepy_core.cpp)

# add target
pybind11_add_module(splinepy_core ${PYSPLINEPY_SRCS})
target_include_directories(splinepy_core PRIVATE "${PROJECT_SOURCE_DIR}/cpp")
target_compile_definitions(splinepy_core PRIVATE ${DEFS})
target_compile_options(splinepy_core
    PRIVATE
        ${SPLINEPY_FLAGS}
        ${SPLINEPY_OPTIMIZATION_FLAGS})

# link splinepy
target_link_libraries(splinepy_core PRIVATE splinepy)

# link splinelib
if(SPLINEPY_COMPILE_SPLINELIB)
    target_link_libraries(splinepy_core PRIVATE splines input_output)
else()
    find_package(SplineLib REQUIRED)
    target_link_libraries(splinepy_core
        PRIVATE
            SplineLib::splines
            SplineLib::input_output
    )
endif()

# link bezman
target_link_libraries(splinepy_core PRIVATE bezman)

# link fitting
target_link_libraries(splinepy_core PRIVATE fitting)

# link napf
target_link_libraries(splinepy_core PRIVATE napf)
