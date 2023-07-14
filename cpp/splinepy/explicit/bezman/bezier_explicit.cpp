#include <bezman/src/bezier_group.hpp>
#include <bezman/src/bezier_spline.hpp>
#include <bezman/src/vertices.hpp>
template class bezman::BezierSpline<1, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<2, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<3, bezman::VertexView<double>, double>;
#ifdef SPLINEPY_MORE
template class bezman::BezierSpline<4, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<5, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<6, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<7, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<8, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<9, bezman::VertexView<double>, double>;
template class bezman::BezierSpline<10, bezman::VertexView<double>, double>;
#endif
