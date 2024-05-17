/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#ifdef SPLINEPY_BUILD_EXPLICIT
extern template class splinepy::splines::RationalBezier<1, 1>;
extern template class splinepy::splines::RationalBezier<1, 2>;
extern template class splinepy::splines::RationalBezier<1, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::RationalBezier<1, 4>;
extern template class splinepy::splines::RationalBezier<1, 5>;
extern template class splinepy::splines::RationalBezier<1, 6>;
extern template class splinepy::splines::RationalBezier<1, 7>;
extern template class splinepy::splines::RationalBezier<1, 8>;
extern template class splinepy::splines::RationalBezier<1, 9>;
extern template class splinepy::splines::RationalBezier<1, 10>;
#endif
extern template class splinepy::splines::RationalBezier<2, 1>;
extern template class splinepy::splines::RationalBezier<2, 2>;
extern template class splinepy::splines::RationalBezier<2, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::RationalBezier<2, 4>;
extern template class splinepy::splines::RationalBezier<2, 5>;
extern template class splinepy::splines::RationalBezier<2, 6>;
extern template class splinepy::splines::RationalBezier<2, 7>;
extern template class splinepy::splines::RationalBezier<2, 8>;
extern template class splinepy::splines::RationalBezier<2, 9>;
extern template class splinepy::splines::RationalBezier<2, 10>;
#endif
extern template class splinepy::splines::RationalBezier<3, 1>;
extern template class splinepy::splines::RationalBezier<3, 2>;
extern template class splinepy::splines::RationalBezier<3, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::RationalBezier<3, 4>;
extern template class splinepy::splines::RationalBezier<3, 5>;
extern template class splinepy::splines::RationalBezier<3, 6>;
extern template class splinepy::splines::RationalBezier<3, 7>;
extern template class splinepy::splines::RationalBezier<3, 8>;
extern template class splinepy::splines::RationalBezier<3, 9>;
extern template class splinepy::splines::RationalBezier<3, 10>;
extern template class splinepy::splines::RationalBezier<4, 1>;
extern template class splinepy::splines::RationalBezier<4, 2>;
extern template class splinepy::splines::RationalBezier<4, 3>;
extern template class splinepy::splines::RationalBezier<4, 4>;
extern template class splinepy::splines::RationalBezier<4, 5>;
extern template class splinepy::splines::RationalBezier<4, 6>;
extern template class splinepy::splines::RationalBezier<4, 7>;
extern template class splinepy::splines::RationalBezier<4, 8>;
extern template class splinepy::splines::RationalBezier<4, 9>;
extern template class splinepy::splines::RationalBezier<4, 10>;

extern template class splinepy::splines::RationalBezier<5, 1>;
extern template class splinepy::splines::RationalBezier<5, 2>;
extern template class splinepy::splines::RationalBezier<5, 3>;
extern template class splinepy::splines::RationalBezier<5, 4>;
extern template class splinepy::splines::RationalBezier<5, 5>;
extern template class splinepy::splines::RationalBezier<5, 6>;
extern template class splinepy::splines::RationalBezier<5, 7>;
extern template class splinepy::splines::RationalBezier<5, 8>;
extern template class splinepy::splines::RationalBezier<5, 9>;
extern template class splinepy::splines::RationalBezier<5, 10>;

extern template class splinepy::splines::RationalBezier<6, 1>;
extern template class splinepy::splines::RationalBezier<6, 2>;
extern template class splinepy::splines::RationalBezier<6, 3>;
extern template class splinepy::splines::RationalBezier<6, 4>;
extern template class splinepy::splines::RationalBezier<6, 5>;
extern template class splinepy::splines::RationalBezier<6, 6>;
extern template class splinepy::splines::RationalBezier<6, 7>;
extern template class splinepy::splines::RationalBezier<6, 8>;
extern template class splinepy::splines::RationalBezier<6, 9>;
extern template class splinepy::splines::RationalBezier<6, 10>;

extern template class splinepy::splines::RationalBezier<7, 1>;
extern template class splinepy::splines::RationalBezier<7, 2>;
extern template class splinepy::splines::RationalBezier<7, 3>;
extern template class splinepy::splines::RationalBezier<7, 4>;
extern template class splinepy::splines::RationalBezier<7, 5>;
extern template class splinepy::splines::RationalBezier<7, 6>;
extern template class splinepy::splines::RationalBezier<7, 7>;
extern template class splinepy::splines::RationalBezier<7, 8>;
extern template class splinepy::splines::RationalBezier<7, 9>;
extern template class splinepy::splines::RationalBezier<7, 10>;

extern template class splinepy::splines::RationalBezier<8, 1>;
extern template class splinepy::splines::RationalBezier<8, 2>;
extern template class splinepy::splines::RationalBezier<8, 3>;
extern template class splinepy::splines::RationalBezier<8, 4>;
extern template class splinepy::splines::RationalBezier<8, 5>;
extern template class splinepy::splines::RationalBezier<8, 6>;
extern template class splinepy::splines::RationalBezier<8, 7>;
extern template class splinepy::splines::RationalBezier<8, 8>;
extern template class splinepy::splines::RationalBezier<8, 9>;
extern template class splinepy::splines::RationalBezier<8, 10>;

extern template class splinepy::splines::RationalBezier<9, 1>;
extern template class splinepy::splines::RationalBezier<9, 2>;
extern template class splinepy::splines::RationalBezier<9, 3>;
extern template class splinepy::splines::RationalBezier<9, 4>;
extern template class splinepy::splines::RationalBezier<9, 5>;
extern template class splinepy::splines::RationalBezier<9, 6>;
extern template class splinepy::splines::RationalBezier<9, 7>;
extern template class splinepy::splines::RationalBezier<9, 8>;
extern template class splinepy::splines::RationalBezier<9, 9>;
extern template class splinepy::splines::RationalBezier<9, 10>;

extern template class splinepy::splines::RationalBezier<10, 1>;
extern template class splinepy::splines::RationalBezier<10, 2>;
extern template class splinepy::splines::RationalBezier<10, 3>;
extern template class splinepy::splines::RationalBezier<10, 4>;
extern template class splinepy::splines::RationalBezier<10, 5>;
extern template class splinepy::splines::RationalBezier<10, 6>;
extern template class splinepy::splines::RationalBezier<10, 7>;
extern template class splinepy::splines::RationalBezier<10, 8>;
extern template class splinepy::splines::RationalBezier<10, 9>;
extern template class splinepy::splines::RationalBezier<10, 10>;
#endif
#endif
