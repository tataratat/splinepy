#pragma once

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class splinepy::splines::BSpline<1>;
extern template class splinepy::splines::BSpline<2>;
extern template class splinepy::splines::BSpline<3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::BSpline<4>;
extern template class splinepy::splines::BSpline<5>;
extern template class splinepy::splines::BSpline<6>;
extern template class splinepy::splines::BSpline<7>;
extern template class splinepy::splines::BSpline<8>;
extern template class splinepy::splines::BSpline<9>;
extern template class splinepy::splines::BSpline<10>;
#endif
#endif
