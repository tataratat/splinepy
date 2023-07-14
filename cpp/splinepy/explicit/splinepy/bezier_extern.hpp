#pragma once
#include <splinepy/splines/bezier.hpp>

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class splinepy::splines::Bezier<1>;
extern template class splinepy::splines::Bezier<2>;
extern template class splinepy::splines::Bezier<3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Bezier<4>;
extern template class splinepy::splines::Bezier<5>;
extern template class splinepy::splines::Bezier<6>;
extern template class splinepy::splines::Bezier<7>;
extern template class splinepy::splines::Bezier<8>;
extern template class splinepy::splines::Bezier<9>;
extern template class splinepy::splines::Bezier<10>;
#endif
#endif
