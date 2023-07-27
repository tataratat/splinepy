#include "splinepy/splines/bezier.hpp"

template class splinepy::splines::Bezier<2, 1>;
template class splinepy::splines::Bezier<2, 2>;
template class splinepy::splines::Bezier<2, 3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::Bezier<2, 4>;
template class splinepy::splines::Bezier<2, 5>;
template class splinepy::splines::Bezier<2, 6>;
template class splinepy::splines::Bezier<2, 7>;
template class splinepy::splines::Bezier<2, 8>;
template class splinepy::splines::Bezier<2, 9>;
template class splinepy::splines::Bezier<2, 10>;
#endif
