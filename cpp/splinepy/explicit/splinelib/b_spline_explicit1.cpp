#include "Sources/Splines/b_spline.hpp"
template class splinelib::sources::splines::BSpline<1, 1>;
template class splinelib::sources::splines::BSpline<1, 2>;
template class splinelib::sources::splines::BSpline<1, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::BSpline<1, 4>;
template class splinelib::sources::splines::BSpline<1, 5>;
template class splinelib::sources::splines::BSpline<1, 6>;
template class splinelib::sources::splines::BSpline<1, 7>;
template class splinelib::sources::splines::BSpline<1, 8>;
template class splinelib::sources::splines::BSpline<1, 9>;
template class splinelib::sources::splines::BSpline<1, 10>;
#endif
