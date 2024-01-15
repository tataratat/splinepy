# reader, export helper, and python module
set(PYSPLINEPY_SRCS
    py_coordinate_pointers.cpp
    py_knot_insertion_matrix.cpp
    py_knot_vector.cpp
    py_multipatch.cpp
    py_spline.cpp
    py_spline_exporter.cpp
    py_spline_extensions.cpp
    py_spline_reader.cpp)

set(PYSPLINEPY_MODULE_SRCS splinepy_core.cpp)

# add target - build an object only with python bindings it is easier to link
# this to an external library, if you prefer to have access to PySpline, instead
# of py::object
add_library(splinepy_python ${SPLINEPY_LIB_TYPE} ${PYSPLINEPY_SRCS})
add_library(splinepy::splinepy_python ALIAS splinepy_python)
target_link_libraries(splinepy_python PUBLIC pybind11::module
                                             splinepy::splinepy pybind11::lto)
target_compile_definitions(splinepy_python PRIVATE ${SPLINEPY_DEFS})
target_compile_options(splinepy_python PRIVATE ${SPLINEPY_FLAGS})

# add module
pybind11_add_module(splinepy_core ${PYSPLINEPY_MODULE_SRCS})
target_compile_definitions(splinepy_core PRIVATE ${SPLINEPY_DEFS})
target_compile_options(splinepy_core PRIVATE ${SPLINEPY_FLAGS})

# link splinepy - all the other dependencies should propagte from splinepy
target_link_libraries(splinepy_core PUBLIC splinepy_python)

# install
install(
  TARGETS splinepy_python
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION "${incl_dest}")

install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE "${namespace}"
  DESTINATION "${cfg_dest}")

install(
  TARGETS splinepy_core
  DESTINATION splinepy
  COMPONENT PythonModule)
