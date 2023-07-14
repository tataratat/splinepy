#pragma once
// related headers
#include <bezman/src/bezier_group.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
#ifdef SPLINEPY_MORE
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::VertexView<double>, double>;
#endif
#endif
