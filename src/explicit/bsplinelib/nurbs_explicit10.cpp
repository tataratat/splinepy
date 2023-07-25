#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<10, 1>;
template class bsplinelib::splines::Nurbs<10, 2>;
template class bsplinelib::splines::Nurbs<10, 3>;
template class bsplinelib::splines::Nurbs<10, 4>;
template class bsplinelib::splines::Nurbs<10, 5>;
template class bsplinelib::splines::Nurbs<10, 6>;
template class bsplinelib::splines::Nurbs<10, 7>;
template class bsplinelib::splines::Nurbs<10, 8>;
template class bsplinelib::splines::Nurbs<10, 9>;
template class bsplinelib::splines::Nurbs<10, 10>;
#endif
