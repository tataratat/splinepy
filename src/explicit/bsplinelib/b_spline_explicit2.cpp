#include <BSplineLib/Splines/b_spline.hpp>

#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_vector_space.hpp"

template class bsplinelib::splines::BSpline<2, 1>;
template class bsplinelib::splines::BSpline<2, 2>;
template class bsplinelib::splines::BSpline<2, 3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::splines::BSpline<2, 4>;
template class bsplinelib::splines::BSpline<2, 5>;
template class bsplinelib::splines::BSpline<2, 6>;
template class bsplinelib::splines::BSpline<2, 7>;
template class bsplinelib::splines::BSpline<2, 8>;
template class bsplinelib::splines::BSpline<2, 9>;
template class bsplinelib::splines::BSpline<2, 10>;
#endif
