#include "BSplineLib/Splines/nurbs.hpp"
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
