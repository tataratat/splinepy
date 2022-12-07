#include "Sources/Splines/nurbs.hpp"
template class splinelib::sources::splines::Nurbs<2, 1>;
template class splinelib::sources::splines::Nurbs<2, 2>;
template class splinelib::sources::splines::Nurbs<2, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::Nurbs<2, 4>;
template class splinelib::sources::splines::Nurbs<2, 5>;
template class splinelib::sources::splines::Nurbs<2, 6>;
template class splinelib::sources::splines::Nurbs<2, 7>;
template class splinelib::sources::splines::Nurbs<2, 8>;
template class splinelib::sources::splines::Nurbs<2, 9>;
template class splinelib::sources::splines::Nurbs<2, 10>;
#endif
