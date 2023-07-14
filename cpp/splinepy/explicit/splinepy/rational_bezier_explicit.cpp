#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/rational_bezier.inl>
template class splinepy::splines::RationalBezier<1>;
template class splinepy::splines::RationalBezier<2>;
template class splinepy::splines::RationalBezier<3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::RationalBezier<4>;
template class splinepy::splines::RationalBezier<5>;
template class splinepy::splines::RationalBezier<6>;
template class splinepy::splines::RationalBezier<7>;
template class splinepy::splines::RationalBezier<8>;
template class splinepy::splines::RationalBezier<9>;
template class splinepy::splines::RationalBezier<10>;
#endif
