#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

template class bsplinelib::splines::Nurbs<3, 1>;
template class bsplinelib::splines::Nurbs<3, 2>;
template class bsplinelib::splines::Nurbs<3, 3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<3, 4>;
template class bsplinelib::splines::Nurbs<3, 5>;
template class bsplinelib::splines::Nurbs<3, 6>;
template class bsplinelib::splines::Nurbs<3, 7>;
template class bsplinelib::splines::Nurbs<3, 8>;
template class bsplinelib::splines::Nurbs<3, 9>;
template class bsplinelib::splines::Nurbs<3, 10>;
#endif
