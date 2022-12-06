#pragma once
// related headers
#include <bezman/src/bezier_group.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class bezman::RationalBezierSpline<1, double, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<3>, double>;
#ifdef SPLINEPY_MORE
extern template class bezman::RationalBezierSpline<1, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<1, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<1, bezman::Point<10>, double>;
#endif
extern template class bezman::RationalBezierSpline<2, double, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<3>, double>;
#ifdef SPLINEPY_MORE
extern template class bezman::RationalBezierSpline<2, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<2, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<2, bezman::Point<10>, double>;
#endif
extern template class bezman::RationalBezierSpline<3, double, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<3>, double>;
#ifdef SPLINEPY_MORE
extern template class bezman::RationalBezierSpline<3, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<3, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<3, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<4, double, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<4, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<4, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<5, double, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<5, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<5, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<6, double, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<6, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<6, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<7, double, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<7, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<7, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<8, double, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<8, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<8, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<9, double, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<2>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<3>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<4>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<5>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<6>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<7>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<8>, double>;
extern template class bezman::RationalBezierSpline<9, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<9, bezman::Point<10>, double>;

extern template class bezman::RationalBezierSpline<10, double, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<2>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<3>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<4>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<5>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<6>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<7>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<8>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<9>, double>;
extern template class bezman::
    RationalBezierSpline<10, bezman::Point<10>, double>;
#endif

#endif
