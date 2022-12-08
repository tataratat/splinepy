#include "Sources/Splines/b_spline.hpp"
template class splinelib::sources::splines::BSpline<3, 1>;
template class splinelib::sources::splines::BSpline<3, 2>;
template class splinelib::sources::splines::BSpline<3, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::BSpline<3, 4>;
template class splinelib::sources::splines::BSpline<3, 5>;
template class splinelib::sources::splines::BSpline<3, 6>;
template class splinelib::sources::splines::BSpline<3, 7>;
template class splinelib::sources::splines::BSpline<3, 8>;
template class splinelib::sources::splines::BSpline<3, 9>;
template class splinelib::sources::splines::BSpline<3, 10>;
#endif
