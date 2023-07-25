#pragma once

extern template class bsplinelib::vector_spaces::WeightedVectorSpace<1>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<2>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<3>;
#ifdef SPLINEPY_MORE
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<4>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<5>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<6>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<7>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<8>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<9>;
extern template class bsplinelib::vector_spaces::WeightedVectorSpace<10>;
#endif
