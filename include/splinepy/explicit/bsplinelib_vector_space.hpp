#pragma once

extern template class bsplinelib::vector_spaces::VectorSpace<1>;
extern template class bsplinelib::vector_spaces::VectorSpace<2>;
extern template class bsplinelib::vector_spaces::VectorSpace<3>;
#ifdef SPLINEPY_MORE
extern template class bsplinelib::vector_spaces::VectorSpace<4>;
extern template class bsplinelib::vector_spaces::VectorSpace<5>;
extern template class bsplinelib::vector_spaces::VectorSpace<6>;
extern template class bsplinelib::vector_spaces::VectorSpace<7>;
extern template class bsplinelib::vector_spaces::VectorSpace<8>;
extern template class bsplinelib::vector_spaces::VectorSpace<9>;
extern template class bsplinelib::vector_spaces::VectorSpace<10>;
#endif
