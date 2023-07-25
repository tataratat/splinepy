#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<6, 1>;
template class bsplinelib::splines::Nurbs<6, 2>;
template class bsplinelib::splines::Nurbs<6, 3>;
template class bsplinelib::splines::Nurbs<6, 4>;
template class bsplinelib::splines::Nurbs<6, 5>;
template class bsplinelib::splines::Nurbs<6, 6>;
template class bsplinelib::splines::Nurbs<6, 7>;
template class bsplinelib::splines::Nurbs<6, 8>;
template class bsplinelib::splines::Nurbs<6, 9>;
template class bsplinelib::splines::Nurbs<6, 10>;
#endif
