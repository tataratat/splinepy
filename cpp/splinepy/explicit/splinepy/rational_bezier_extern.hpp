#pragma once
#include <splinepy/splines/rational_bezier.hpp>

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class splinepy::splines::RationalBezier<1>;
extern template class splinepy::splines::RationalBezier<2>;
extern template class splinepy::splines::RationalBezier<3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::RationalBezier<4>;
extern template class splinepy::splines::RationalBezier<5>;
extern template class splinepy::splines::RationalBezier<6>;
extern template class splinepy::splines::RationalBezier<7>;
extern template class splinepy::splines::RationalBezier<8>;
extern template class splinepy::splines::RationalBezier<9>;
extern template class splinepy::splines::RationalBezier<10>;
#endif
#endif
