#include "Sources/Splines/nurbs.hpp"
template class splinelib::sources::splines::Nurbs<3, 1>;
template class splinelib::sources::splines::Nurbs<3, 2>;
template class splinelib::sources::splines::Nurbs<3, 3>;
#ifdef SPLINEPY_MORE
template class splinelib::sources::splines::Nurbs<3, 4>;
template class splinelib::sources::splines::Nurbs<3, 5>;
template class splinelib::sources::splines::Nurbs<3, 6>;
template class splinelib::sources::splines::Nurbs<3, 7>;
template class splinelib::sources::splines::Nurbs<3, 8>;
template class splinelib::sources::splines::Nurbs<3, 9>;
template class splinelib::sources::splines::Nurbs<3, 10>;
#endif
