#pragma once
#include "BSplineLib/Splines/nurbs.hpp"

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class bsplinelib::splines::Nurbs<1>;
extern template class bsplinelib::splines::Nurbs<2>;
extern template class bsplinelib::splines::Nurbs<3>;
#ifdef SPLINEPY_MORE
extern template class bsplinelib::splines::Nurbs<4>;
extern template class bsplinelib::splines::Nurbs<5>;
extern template class bsplinelib::splines::Nurbs<6>;
extern template class bsplinelib::splines::Nurbs<7>;
extern template class bsplinelib::splines::Nurbs<8>;
extern template class bsplinelib::splines::Nurbs<9>;
extern template class bsplinelib::splines::Nurbs<10>;
#endif

#endif
