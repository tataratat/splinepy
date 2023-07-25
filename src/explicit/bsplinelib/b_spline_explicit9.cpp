#include <BSplineLib/Splines/b_spline.hpp>

#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::BSpline<9, 1>;
template class bsplinelib::splines::BSpline<9, 2>;
template class bsplinelib::splines::BSpline<9, 3>;
template class bsplinelib::splines::BSpline<9, 4>;
template class bsplinelib::splines::BSpline<9, 5>;
template class bsplinelib::splines::BSpline<9, 6>;
template class bsplinelib::splines::BSpline<9, 7>;
template class bsplinelib::splines::BSpline<9, 8>;
template class bsplinelib::splines::BSpline<9, 9>;
template class bsplinelib::splines::BSpline<9, 10>;
#endif
