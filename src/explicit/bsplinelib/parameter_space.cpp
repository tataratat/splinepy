#include <BSplineLib/ParameterSpaces/parameter_space.hpp>

template class bsplinelib::parameter_spaces::ParameterSpace<1>;
template class bsplinelib::parameter_spaces::ParameterSpace<2>;
template class bsplinelib::parameter_spaces::ParameterSpace<3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::parameter_spaces::ParameterSpace<4>;
template class bsplinelib::parameter_spaces::ParameterSpace<5>;
template class bsplinelib::parameter_spaces::ParameterSpace<6>;
template class bsplinelib::parameter_spaces::ParameterSpace<7>;
template class bsplinelib::parameter_spaces::ParameterSpace<8>;
template class bsplinelib::parameter_spaces::ParameterSpace<9>;
template class bsplinelib::parameter_spaces::ParameterSpace<10>;
#endif
