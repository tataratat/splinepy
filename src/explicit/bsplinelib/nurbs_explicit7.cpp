#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<7, 1>;
template class bsplinelib::splines::Nurbs<7, 2>;
template class bsplinelib::splines::Nurbs<7, 3>;
template class bsplinelib::splines::Nurbs<7, 4>;
template class bsplinelib::splines::Nurbs<7, 5>;
template class bsplinelib::splines::Nurbs<7, 6>;
template class bsplinelib::splines::Nurbs<7, 7>;
template class bsplinelib::splines::Nurbs<7, 8>;
template class bsplinelib::splines::Nurbs<7, 9>;
template class bsplinelib::splines::Nurbs<7, 10>;
#endif
