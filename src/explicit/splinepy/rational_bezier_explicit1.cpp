#include "splinepy/splines/rational_bezier.hpp"

template class splinepy::splines::RationalBezier<1, 1>;
template class splinepy::splines::RationalBezier<1, 2>;
template class splinepy::splines::RationalBezier<1, 3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::RationalBezier<1, 4>;
template class splinepy::splines::RationalBezier<1, 5>;
template class splinepy::splines::RationalBezier<1, 6>;
template class splinepy::splines::RationalBezier<1, 7>;
template class splinepy::splines::RationalBezier<1, 8>;
template class splinepy::splines::RationalBezier<1, 9>;
template class splinepy::splines::RationalBezier<1, 10>;
#endif
