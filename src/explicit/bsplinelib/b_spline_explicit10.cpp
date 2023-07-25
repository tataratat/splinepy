#include <BSplineLib/Splines/b_spline.hpp>

#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::BSpline<10, 1>;
template class bsplinelib::splines::BSpline<10, 2>;
template class bsplinelib::splines::BSpline<10, 3>;
template class bsplinelib::splines::BSpline<10, 4>;
template class bsplinelib::splines::BSpline<10, 5>;
template class bsplinelib::splines::BSpline<10, 6>;
template class bsplinelib::splines::BSpline<10, 7>;
template class bsplinelib::splines::BSpline<10, 8>;
template class bsplinelib::splines::BSpline<10, 9>;
template class bsplinelib::splines::BSpline<10, 10>;
#endif
