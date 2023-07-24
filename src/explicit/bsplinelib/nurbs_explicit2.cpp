#include "BSplineLib/Splines/nurbs.hpp"
template class bsplinelib::splines::Nurbs<2, 1>;
template class bsplinelib::splines::Nurbs<2, 2>;
template class bsplinelib::splines::Nurbs<2, 3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<2, 4>;
template class bsplinelib::splines::Nurbs<2, 5>;
template class bsplinelib::splines::Nurbs<2, 6>;
template class bsplinelib::splines::Nurbs<2, 7>;
template class bsplinelib::splines::Nurbs<2, 8>;
template class bsplinelib::splines::Nurbs<2, 9>;
template class bsplinelib::splines::Nurbs<2, 10>;
#endif
