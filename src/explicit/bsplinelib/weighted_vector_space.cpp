#include <BSplineLib/VectorSpaces/weighted_vector_space.hpp>

template class bsplinelib::vector_spaces::WeightedVectorSpace<1>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<2>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<3>;
#ifdef SPLINEPY_MORE
template class bsplinelib::vector_spaces::WeightedVectorSpace<4>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<5>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<6>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<7>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<8>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<9>;
template class bsplinelib::vector_spaces::WeightedVectorSpace<10>;
#endif
