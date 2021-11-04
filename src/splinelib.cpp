#include <spline_reader.hpp>
#include <bspline.hpp>
#include <nurbs.hpp>

PYBIND11_MODULE(splinelibpy, m) {

    py::class_<PyBSpline<1, 2>>(m, "BSplineCurve2D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyBSpline<1, 2>::p_knot_vectors)
        .def_readwrite("degrees", &PyBSpline<1, 2>::p_degrees)
        .def_readwrite("control_points", &PyBSpline<1, 2>::p_control_points)
        .def("evaluate", &PyBSpline<1, 2>::evaluate, py::arg("queries"))
        .def("derivative", &PyBSpline<1, 2>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyBSpline<1, 2>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyBSpline<1, 2>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyBSpline<1, 2>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyBSpline<1, 2>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyBSpline<1, 2>::sample, py::arg("resoultion"))
        .def("fit_curve", &PyBSpline<1, 2>::fit_curve, py::arg("points"),
                                                       py::arg("degree"),
                                                       py::arg("num_control_points"),
                                                       py::arg("centripetal"),
                                                       py::arg("knot_vector"))
        .def("interpolate_curve", &PyBSpline<1, 2>::interpolate_curve, py::arg("points"),
                                                                       py::arg("degree"),
                                                                       py::arg("centripetal"))
        .def("approximate_curve", &PyBSpline<1, 2>::approximate_curve, py::arg("points"),
                                                                       py::arg("degree"),
                                                                       py::arg("num_control_points"),
                                                                       py::arg("centripetal"))
        .def("write_iges", &PyBSpline<1, 2>::write_iges, py::arg("fname"))
        .def("write_xml", &PyBSpline<1, 2>::write_xml, py::arg("fname"))
        .def("write_irit", &PyBSpline<1, 2>::write_irit, py::arg("fname"))
        .def("update_c_", &PyBSpline<1, 2>::update_c)
        .def("update_p_", &PyBSpline<1, 2>::update_p)
        ;

    py::class_<PyBSpline<1, 3>>(m, "BSplineCurve3D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyBSpline<1, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyBSpline<1, 3>::p_degrees)
        .def_readwrite("control_points", &PyBSpline<1, 3>::p_control_points)
        .def("evaluate", &PyBSpline<1, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyBSpline<1, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyBSpline<1, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyBSpline<1, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyBSpline<1, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyBSpline<1, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyBSpline<1, 3>::sample, py::arg("resoultion"))
        .def("fit_curve", &PyBSpline<1, 3>::fit_curve, py::arg("points"),
                                                       py::arg("degree"),
                                                       py::arg("num_control_points"),
                                                       py::arg("centripetal"),
                                                       py::arg("knot_vector"))
        .def("interpolate_curve", &PyBSpline<1, 3>::interpolate_curve, py::arg("points"),
                                                                       py::arg("degree"),
                                                                       py::arg("centripetal"))
        .def("approximate_curve", &PyBSpline<1, 3>::approximate_curve, py::arg("points"),
                                                                       py::arg("degree"),
                                                                       py::arg("num_control_points"),
                                                                       py::arg("centripetal"))
        .def("write_iges", &PyBSpline<1, 3>::write_iges, py::arg("fname"))
        .def("write_xml", &PyBSpline<1, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyBSpline<1, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyBSpline<1, 3>::update_c)
        .def("update_p_", &PyBSpline<1, 3>::update_p)
        ;

    py::class_<PyBSpline<2, 2>>(m, "BSplineSurface2D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyBSpline<2, 2>::p_knot_vectors)
        .def_readwrite("degrees", &PyBSpline<2, 2>::p_degrees)
        .def_readwrite("control_points", &PyBSpline<2, 2>::p_control_points)
        .def("evaluate", &PyBSpline<2, 2>::evaluate, py::arg("queries"))
        .def("derivative", &PyBSpline<2, 2>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyBSpline<2, 2>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyBSpline<2, 2>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyBSpline<2, 2>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyBSpline<2, 2>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyBSpline<2, 2>::sample, py::arg("resoultion"))
        .def("fit_surface", &PyBSpline<2, 2>::fit_surface, py::arg("points"),
                                                           py::arg("size_u"),
                                                           py::arg("size_v"),
                                                           py::arg("degree_u"),
                                                           py::arg("degree_v"),
                                                           py::arg("centripetal"))
        .def("interpolate_surface", &PyBSpline<2, 2>::interpolate_surface, py::arg("points"),
                                                                           py::arg("size_u"),
                                                                           py::arg("size_v"),
                                                                           py::arg("degree_u"),
                                                                           py::arg("degree_v"),
                                                                           py::arg("centripetal"))
        .def("write_iges", &PyBSpline<2, 2>::write_iges, py::arg("fname"))
        .def("write_xml", &PyBSpline<2, 2>::write_xml, py::arg("fname"))
        .def("write_irit", &PyBSpline<2, 2>::write_irit, py::arg("fname"))
        .def("update_c_", &PyBSpline<2, 2>::update_c)
        .def("update_p_", &PyBSpline<2, 2>::update_p)
        ;

    py::class_<PyBSpline<2, 3>>(m, "BSplineSurface3D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyBSpline<2, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyBSpline<2, 3>::p_degrees)
        .def_readwrite("control_points", &PyBSpline<2, 3>::p_control_points)
        .def("evaluate", &PyBSpline<2, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyBSpline<2, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyBSpline<2, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyBSpline<2, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyBSpline<2, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyBSpline<2, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))  
        .def("sample", &PyBSpline<2, 3>::sample, py::arg("resoultion"))
        .def("fit_surface", &PyBSpline<2, 3>::fit_surface, py::arg("points"),
                                                           py::arg("size_u"),
                                                           py::arg("size_v"),
                                                           py::arg("degree_u"),
                                                           py::arg("degree_v"),
                                                           py::arg("centripetal"))
        .def("interpolate_surface", &PyBSpline<2, 3>::interpolate_surface, py::arg("points"),
                                                                           py::arg("size_u"),
                                                                           py::arg("size_v"),
                                                                           py::arg("degree_u"),
                                                                           py::arg("degree_v"),
                                                                           py::arg("centripetal"))
        .def("write_iges", &PyBSpline<2, 3>::write_iges, py::arg("fname"))
        .def("write_xml", &PyBSpline<2, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyBSpline<2, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyBSpline<2, 3>::update_c)
        .def("update_p_", &PyBSpline<2, 3>::update_p)
        ;

    py::class_<PyBSpline<3, 3>>(m, "BSplineSolid")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyBSpline<3, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyBSpline<3, 3>::p_degrees)
        .def_readwrite("control_points", &PyBSpline<3, 3>::p_control_points)
        .def("evaluate", &PyBSpline<3, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyBSpline<3, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyBSpline<3, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyBSpline<3, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyBSpline<3, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyBSpline<3, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyBSpline<3, 3>::sample, py::arg("resoultion"))
        .def("write_xml", &PyBSpline<3, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyBSpline<3, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyBSpline<3, 3>::update_c)
        .def("update_p_", &PyBSpline<3, 3>::update_p)
        ;

    py::class_<PyNurbs<1, 2>>(m, "NurbsCurve2D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyNurbs<1, 2>::p_knot_vectors)
        .def_readwrite("degrees", &PyNurbs<1, 2>::p_degrees)
        .def_readwrite("control_points", &PyNurbs<1, 2>::p_control_points)
        .def_readwrite("weights", &PyNurbs<1, 2>::p_weights)
        .def("evaluate", &PyNurbs<1, 2>::evaluate, py::arg("queries"))
        .def("derivative", &PyNurbs<1, 2>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyNurbs<1, 2>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyNurbs<1, 2>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyNurbs<1, 2>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyNurbs<1, 2>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyNurbs<1, 2>::sample, py::arg("resoultion"))
        .def("write_iges", &PyNurbs<1, 2>::write_iges, py::arg("fname"))
        .def("write_xml", &PyNurbs<1, 2>::write_xml, py::arg("fname"))
        .def("write_irit", &PyNurbs<1, 2>::write_irit, py::arg("fname"))
        .def("update_c_", &PyNurbs<1, 2>::update_c)
        .def("update_p_", &PyNurbs<1, 2>::update_p)
        ;

    py::class_<PyNurbs<1, 3>>(m, "NurbsCurve3D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyNurbs<1, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyNurbs<1, 3>::p_degrees)
        .def_readwrite("control_points", &PyNurbs<1, 3>::p_control_points)
        .def_readwrite("weights", &PyNurbs<1, 3>::p_weights)
        .def("evaluate", &PyNurbs<1, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyNurbs<1, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyNurbs<1, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyNurbs<1, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyNurbs<1, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyNurbs<1, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyNurbs<1, 3>::sample, py::arg("resoultion"))
        .def("write_iges", &PyNurbs<1, 3>::write_iges, py::arg("fname"))
        .def("write_xml", &PyNurbs<1, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyNurbs<1, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyNurbs<1, 3>::update_c)
        .def("update_p_", &PyNurbs<1, 3>::update_p)
        ;

    py::class_<PyNurbs<2, 2>>(m, "NurbsSurface2D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyNurbs<2, 2>::p_knot_vectors)
        .def_readwrite("degrees", &PyNurbs<2, 2>::p_degrees)
        .def_readwrite("control_points", &PyNurbs<2, 2>::p_control_points)
        .def_readwrite("weights", &PyNurbs<2, 2>::p_weights)
        .def("evaluate", &PyNurbs<2, 2>::evaluate, py::arg("queries"))
        .def("derivative", &PyNurbs<2, 2>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyNurbs<2, 2>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyNurbs<2, 2>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyNurbs<2, 2>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyNurbs<2, 2>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyNurbs<2, 2>::sample, py::arg("resoultion"))
        .def("write_iges", &PyNurbs<2, 2>::write_iges, py::arg("fname"))
        .def("write_xml", &PyNurbs<2, 2>::write_xml, py::arg("fname"))
        .def("write_irit", &PyNurbs<2, 2>::write_irit, py::arg("fname"))
        .def("update_c_", &PyNurbs<2, 2>::update_c)
        .def("update_p_", &PyNurbs<2, 2>::update_p)
        ;

    py::class_<PyNurbs<2, 3>>(m, "NurbsSurface3D")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyNurbs<2, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyNurbs<2, 3>::p_degrees)
        .def_readwrite("control_points", &PyNurbs<2, 3>::p_control_points)
        .def_readwrite("weights", &PyNurbs<2, 3>::p_weights)
        .def("evaluate", &PyNurbs<2, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyNurbs<2, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyNurbs<2, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyNurbs<2, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyNurbs<2, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyNurbs<2, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))  
        .def("sample", &PyNurbs<2, 3>::sample, py::arg("resoultion"))
        .def("write_iges", &PyNurbs<2, 3>::write_iges, py::arg("fname"))
        .def("write_xml", &PyNurbs<2, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyNurbs<2, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyNurbs<2, 3>::update_c)
        .def("update_p_", &PyNurbs<2, 3>::update_p)
        ;

    py::class_<PyNurbs<3, 3>>(m, "NurbsSolid")
        .def(py::init<>())
        .def_readwrite("knot_vectors", &PyNurbs<3, 3>::p_knot_vectors)
        .def_readwrite("degrees", &PyNurbs<3, 3>::p_degrees)
        .def_readwrite("control_points", &PyNurbs<3, 3>::p_control_points)
        .def_readwrite("weights", &PyNurbs<3, 3>::p_weights)
        .def("evaluate", &PyNurbs<3, 3>::evaluate, py::arg("queries"))
        .def("derivative", &PyNurbs<3, 3>::derivative, py::arg("queries"), py::arg("orders"))
        .def("insert_knots", &PyNurbs<3, 3>::insert_knots, py::arg("p_dim"), py::arg("knots"))
        .def("remove_knots", &PyNurbs<3, 3>::remove_knots, py::arg("p_dim"), py::arg("knots"), py::arg("tolerance"))
        .def("elevate_degree", &PyNurbs<3, 3>::elevate_degree, py::arg("p_dim"))
        .def("reduce_degree", &PyNurbs<3, 3>::reduce_degree, py::arg("p_dim"), py::arg("tolerance"))
        .def("sample", &PyNurbs<3, 3>::sample, py::arg("resoultion"))
        .def("write_xml", &PyNurbs<3, 3>::write_xml, py::arg("fname"))
        .def("write_irit", &PyNurbs<3, 3>::write_irit, py::arg("fname"))
        .def("update_c_", &PyNurbs<3, 3>::update_c)
        .def("update_p_", &PyNurbs<3, 3>::update_p)
        ;

    py::class_<SplineReader>(m, "Reader")
        .def(py::init<>())
        .def("read_iges", &SplineReader::read_iges, py::arg("fname"))
        .def("read_xml", &SplineReader::read_xml, py::arg("fname"))
        .def("read_irit", &SplineReader::read_irit, py::arg("fname"))
    ;


}
