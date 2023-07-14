#include <bezman/src/bezier_group.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>
template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<2, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<3, bezman::VertexView<double>, double>;
#ifdef SPLINEPY_MORE
template class bezman::
    RationalBezierSpline<4, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<5, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<6, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<7, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<8, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<9, bezman::VertexView<double>, double>;
template class bezman::
    RationalBezierSpline<10, bezman::VertexView<double>, double>;
#endif
