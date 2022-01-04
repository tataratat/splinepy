#include <spline_reader.hpp>

void init_reader(py::module_ &m) {
    py::class_<SplineReader>(m, "Reader")
        .def(py::init<>())
        .def("read_iges", &SplineReader::read_iges, py::arg("fname"))
        .def("read_xml", &SplineReader::read_xml, py::arg("fname"))
        .def("read_irit", &SplineReader::read_irit, py::arg("fname"))
        ;
}

