#include "BSplineLib/Splines/nurbs.hpp"
template class bsplinelib::splines::Nurbs<1>;
template class bsplinelib::splines::Nurbs<2>;
template class bsplinelib::splines::Nurbs<3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::splines::Nurbs<4>;
template class bsplinelib::splines::Nurbs<5>;
template class bsplinelib::splines::Nurbs<6>;
template class bsplinelib::splines::Nurbs<7>;
template class bsplinelib::splines::Nurbs<8>;
template class bsplinelib::splines::Nurbs<9>;
template class bsplinelib::splines::Nurbs<10>;
#endif
