#include <BSplineLib/VectorSpaces/vector_space.hpp>

template class bsplinelib::vector_spaces::VectorSpace<1>;
template class bsplinelib::vector_spaces::VectorSpace<2>;
template class bsplinelib::vector_spaces::VectorSpace<3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::vector_spaces::VectorSpace<4>;
template class bsplinelib::vector_spaces::VectorSpace<5>;
template class bsplinelib::vector_spaces::VectorSpace<6>;
template class bsplinelib::vector_spaces::VectorSpace<7>;
template class bsplinelib::vector_spaces::VectorSpace<8>;
template class bsplinelib::vector_spaces::VectorSpace<9>;
template class bsplinelib::vector_spaces::VectorSpace<10>;
#endif
