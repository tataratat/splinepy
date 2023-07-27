#include "splinepy/splines/bezier.hpp"

template class splinepy::splines::Bezier<3, 1>;
template class splinepy::splines::Bezier<3, 2>;
template class splinepy::splines::Bezier<3, 3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::Bezier<3, 4>;
template class splinepy::splines::Bezier<3, 5>;
template class splinepy::splines::Bezier<3, 6>;
template class splinepy::splines::Bezier<3, 7>;
template class splinepy::splines::Bezier<3, 8>;
template class splinepy::splines::Bezier<3, 9>;
template class splinepy::splines::Bezier<3, 10>;
#endif
