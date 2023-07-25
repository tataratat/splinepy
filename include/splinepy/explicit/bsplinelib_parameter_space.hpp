#pragma once

extern template class bsplinelib::parameter_spaces::ParameterSpace<1>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<2>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<3>;
#ifdef SPLINEPY_MORE
extern template class bsplinelib::parameter_spaces::ParameterSpace<4>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<5>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<6>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<7>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<8>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<9>;
extern template class bsplinelib::parameter_spaces::ParameterSpace<10>;
#endif
