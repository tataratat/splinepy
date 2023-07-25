#include "splinepy/splines/bezier.hpp"

template class splinepy::splines::Bezier<1, 1>;
template class splinepy::splines::Bezier<1, 2>;
template class splinepy::splines::Bezier<1, 3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::Bezier<1, 4>;
template class splinepy::splines::Bezier<1, 5>;
template class splinepy::splines::Bezier<1, 6>;
template class splinepy::splines::Bezier<1, 7>;
template class splinepy::splines::Bezier<1, 8>;
template class splinepy::splines::Bezier<1, 9>;
template class splinepy::splines::Bezier<1, 10>;
#endif
