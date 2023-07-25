#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<8, 1>;
template class bsplinelib::splines::Nurbs<8, 2>;
template class bsplinelib::splines::Nurbs<8, 3>;
template class bsplinelib::splines::Nurbs<8, 4>;
template class bsplinelib::splines::Nurbs<8, 5>;
template class bsplinelib::splines::Nurbs<8, 6>;
template class bsplinelib::splines::Nurbs<8, 7>;
template class bsplinelib::splines::Nurbs<8, 8>;
template class bsplinelib::splines::Nurbs<8, 9>;
template class bsplinelib::splines::Nurbs<8, 10>;
#endif
