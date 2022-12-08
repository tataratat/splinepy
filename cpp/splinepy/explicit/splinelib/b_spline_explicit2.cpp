#include "Sources/Splines/b_spline.hpp"
template class splinelib::sources::splines::BSpline<2, 1>;
template class splinelib::sources::splines::BSpline<2, 2>;
template class splinelib::sources::splines::BSpline<2, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::BSpline<2, 4>;
template class splinelib::sources::splines::BSpline<2, 5>;
template class splinelib::sources::splines::BSpline<2, 6>;
template class splinelib::sources::splines::BSpline<2, 7>;
template class splinelib::sources::splines::BSpline<2, 8>;
template class splinelib::sources::splines::BSpline<2, 9>;
template class splinelib::sources::splines::BSpline<2, 10>;
#endif
