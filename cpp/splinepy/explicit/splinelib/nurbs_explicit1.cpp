#include "Sources/Splines/nurbs.hpp"
template class splinelib::sources::splines::Nurbs<1, 1>;
template class splinelib::sources::splines::Nurbs<1, 2>;
template class splinelib::sources::splines::Nurbs<1, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::Nurbs<1, 4>;
template class splinelib::sources::splines::Nurbs<1, 5>;
template class splinelib::sources::splines::Nurbs<1, 6>;
template class splinelib::sources::splines::Nurbs<1, 7>;
template class splinelib::sources::splines::Nurbs<1, 8>;
template class splinelib::sources::splines::Nurbs<1, 9>;
template class splinelib::sources::splines::Nurbs<1, 10>;
#endif
