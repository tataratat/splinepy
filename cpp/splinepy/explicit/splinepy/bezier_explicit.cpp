#include <splinepy/splines/bezier.hpp>
#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/bezier.inl>

template class splinepy::splines::Bezier<1>;
template class splinepy::splines::Bezier<2>;
template class splinepy::splines::Bezier<3>;
#ifdef SPLINEPY_MORE
template class splinepy::splines::Bezier<4>;
template class splinepy::splines::Bezier<5>;
template class splinepy::splines::Bezier<6>;
template class splinepy::splines::Bezier<7>;
template class splinepy::splines::Bezier<8>;
template class splinepy::splines::Bezier<9>;
template class splinepy::splines::Bezier<10>;
#endif
