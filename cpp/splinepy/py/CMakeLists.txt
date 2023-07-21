# reader, export helper, and python module
set(PYSPLINEPY_SRCS
        init/coordinate_pointers.cpp
        init/core_spline.cpp
        init/exporter.cpp
        init/fitting.cpp
        init/knot_insertion_matrix.cpp
        init/knot_vector.cpp
        init/multi_patch.cpp
        init/reader.cpp
        init/spline_extensions.cpp
        splinepy_core.cpp)

# add target
pybind11_add_module(splinepy_core ${PYSPLINEPY_SRCS})
target_include_directories(splinepy_core PRIVATE "${PROJECT_SOURCE_DIR}/cpp")
target_compile_definitions(splinepy_core PRIVATE ${SPLINEPY_DEFS})
target_compile_options(splinepy_core
    PRIVATE
        ${SPLINEPY_FLAGS}
        ${SPLINEPY_OPTIMIZATION_FLAGS})

# link splinepy - all the other dependencies should propagte from splinepy
target_link_libraries(splinepy_core PRIVATE splinepy::splinepy)
