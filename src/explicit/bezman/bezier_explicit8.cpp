#include <bezman/src/bezier_spline.hpp>
#include <bezman/src/point.hpp>
#ifdef SPLINEPY_MORE
template class bezman::BezierSpline<8, double, double>;
template class bezman::BezierSpline<8, bezman::Point<2>, double>;
template class bezman::BezierSpline<8, bezman::Point<3>, double>;
template class bezman::BezierSpline<8, bezman::Point<4>, double>;
template class bezman::BezierSpline<8, bezman::Point<5>, double>;
template class bezman::BezierSpline<8, bezman::Point<6>, double>;
template class bezman::BezierSpline<8, bezman::Point<7>, double>;
template class bezman::BezierSpline<8, bezman::Point<8>, double>;
template class bezman::BezierSpline<8, bezman::Point<9>, double>;
template class bezman::BezierSpline<8, bezman::Point<10>, double>;
#endif
