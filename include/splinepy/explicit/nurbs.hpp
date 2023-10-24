#pragma once

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class splinepy::splines::Nurbs<1>;
extern template class splinepy::splines::Nurbs<2>;
extern template class splinepy::splines::Nurbs<3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Nurbs<4>;
extern template class splinepy::splines::Nurbs<5>;
extern template class splinepy::splines::Nurbs<6>;
extern template class splinepy::splines::Nurbs<7>;
extern template class splinepy::splines::Nurbs<8>;
extern template class splinepy::splines::Nurbs<9>;
extern template class splinepy::splines::Nurbs<10>;
#endif
#endif
