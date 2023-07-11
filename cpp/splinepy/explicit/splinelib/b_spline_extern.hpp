#pragma once
#include "BSplineLib/Splines/b_spline.hpp"

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class bsplinelib::splines::BSpline<1>;
extern template class bsplinelib::splines::BSpline<2>;
extern template class bsplinelib::splines::BSpline<3>;
#ifdef SPLINEPY_MORE
extern template class bsplinelib::splines::BSpline<4>;
extern template class bsplinelib::splines::BSpline<5>;
extern template class bsplinelib::splines::BSpline<6>;
extern template class bsplinelib::splines::BSpline<7>;
extern template class bsplinelib::splines::BSpline<8>;
extern template class bsplinelib::splines::BSpline<9>;
extern template class bsplinelib::splines::BSpline<10>;
#endif

#endif
