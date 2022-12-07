#include <bezman/src/bezier_group.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>
template class bezman::RationalBezierSpline<2, double, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<2>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<3>, double>;
#ifdef SPLINEPY_MORE
template class bezman::RationalBezierSpline<2, bezman::Point<4>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<5>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<6>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<7>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<8>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<9>, double>;
template class bezman::RationalBezierSpline<2, bezman::Point<10>, double>;
#endif
