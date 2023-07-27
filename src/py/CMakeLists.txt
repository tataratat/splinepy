# reader, export helper, and python module
set(PYSPLINEPY_SRCS
    py_coordinate_pointers.cpp
    py_fitting.cpp
    py_knot_insertion_matrix.cpp
    py_knot_vector.cpp
    py_multipatch.cpp
    py_spline.cpp
    py_spline_exporter.cpp
    py_spline_extensions.cpp
    py_spline_reader.cpp
    splinepy_core.cpp)

# add target
pybind11_add_module(splinepy_core ${PYSPLINEPY_SRCS})
target_include_directories(splinepy_core
                           PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_compile_definitions(splinepy_core PRIVATE ${SPLINEPY_DEFS})
target_compile_options(splinepy_core PRIVATE ${SPLINEPY_FLAGS}
                                             ${SPLINEPY_OPTIMIZATION_FLAGS})

# link splinepy - all the other dependencies should propagte from splinepy
target_link_libraries(splinepy_core PRIVATE splinepy::splinepy)
