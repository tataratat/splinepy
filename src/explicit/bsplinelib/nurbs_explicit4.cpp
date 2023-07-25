#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<4, 1>;
template class bsplinelib::splines::Nurbs<4, 2>;
template class bsplinelib::splines::Nurbs<4, 3>;
template class bsplinelib::splines::Nurbs<4, 4>;
template class bsplinelib::splines::Nurbs<4, 5>;
template class bsplinelib::splines::Nurbs<4, 6>;
template class bsplinelib::splines::Nurbs<4, 7>;
template class bsplinelib::splines::Nurbs<4, 8>;
template class bsplinelib::splines::Nurbs<4, 9>;
template class bsplinelib::splines::Nurbs<4, 10>;
#endif
