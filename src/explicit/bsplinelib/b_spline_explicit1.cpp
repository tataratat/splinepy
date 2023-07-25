#include <BSplineLib/Splines/b_spline.hpp>

#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_vector_space.hpp"

template class bsplinelib::splines::BSpline<1, 1>;
template class bsplinelib::splines::BSpline<1, 2>;
template class bsplinelib::splines::BSpline<1, 3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::splines::BSpline<1, 4>;
template class bsplinelib::splines::BSpline<1, 5>;
template class bsplinelib::splines::BSpline<1, 6>;
template class bsplinelib::splines::BSpline<1, 7>;
template class bsplinelib::splines::BSpline<1, 8>;
template class bsplinelib::splines::BSpline<1, 9>;
template class bsplinelib::splines::BSpline<1, 10>;
#endif
