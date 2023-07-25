#include <BSplineLib/Splines/nurbs.hpp>

#include "splinepy/explicit/bsplinelib_bspline.hpp"
#include "splinepy/explicit/bsplinelib_parameter_space.hpp"
#include "splinepy/explicit/bsplinelib_weighted_vector_space.hpp"

#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<9, 1>;
template class bsplinelib::splines::Nurbs<9, 2>;
template class bsplinelib::splines::Nurbs<9, 3>;
template class bsplinelib::splines::Nurbs<9, 4>;
template class bsplinelib::splines::Nurbs<9, 5>;
template class bsplinelib::splines::Nurbs<9, 6>;
template class bsplinelib::splines::Nurbs<9, 7>;
template class bsplinelib::splines::Nurbs<9, 8>;
template class bsplinelib::splines::Nurbs<9, 9>;
template class bsplinelib::splines::Nurbs<9, 10>;
#endif
