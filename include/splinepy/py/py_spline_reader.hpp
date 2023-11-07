#pragma once

#include <string>

// pybind11
#include <pybind11/pybind11.h>

namespace splinepy::py {

namespace py = pybind11;

py::list ReadIges(const std::string fname);

} // namespace splinepy::py
