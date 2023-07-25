#include "splinepy/splines/rational_bezier.hpp"

template class splinepy::splines::RationalBezier<3, 1>;
template class splinepy::splines::RationalBezier<3, 2>;
template class splinepy::splines::RationalBezier<3, 3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::RationalBezier<3, 4>;
template class splinepy::splines::RationalBezier<3, 5>;
template class splinepy::splines::RationalBezier<3, 6>;
template class splinepy::splines::RationalBezier<3, 7>;
template class splinepy::splines::RationalBezier<3, 8>;
template class splinepy::splines::RationalBezier<3, 9>;
template class splinepy::splines::RationalBezier<3, 10>;
#endif
