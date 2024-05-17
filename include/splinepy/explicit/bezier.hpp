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
extern template class splinepy::splines::Bezier<1, 1>;
extern template class splinepy::splines::Bezier<1, 2>;
extern template class splinepy::splines::Bezier<1, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Bezier<1, 4>;
extern template class splinepy::splines::Bezier<1, 5>;
extern template class splinepy::splines::Bezier<1, 6>;
extern template class splinepy::splines::Bezier<1, 7>;
extern template class splinepy::splines::Bezier<1, 8>;
extern template class splinepy::splines::Bezier<1, 9>;
extern template class splinepy::splines::Bezier<1, 10>;
#endif
extern template class splinepy::splines::Bezier<2, 1>;
extern template class splinepy::splines::Bezier<2, 2>;
extern template class splinepy::splines::Bezier<2, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Bezier<2, 4>;
extern template class splinepy::splines::Bezier<2, 5>;
extern template class splinepy::splines::Bezier<2, 6>;
extern template class splinepy::splines::Bezier<2, 7>;
extern template class splinepy::splines::Bezier<2, 8>;
extern template class splinepy::splines::Bezier<2, 9>;
extern template class splinepy::splines::Bezier<2, 10>;
#endif
extern template class splinepy::splines::Bezier<3, 1>;
extern template class splinepy::splines::Bezier<3, 2>;
extern template class splinepy::splines::Bezier<3, 3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Bezier<3, 4>;
extern template class splinepy::splines::Bezier<3, 5>;
extern template class splinepy::splines::Bezier<3, 6>;
extern template class splinepy::splines::Bezier<3, 7>;
extern template class splinepy::splines::Bezier<3, 8>;
extern template class splinepy::splines::Bezier<3, 9>;
extern template class splinepy::splines::Bezier<3, 10>;
extern template class splinepy::splines::Bezier<4, 1>;
extern template class splinepy::splines::Bezier<4, 2>;
extern template class splinepy::splines::Bezier<4, 3>;
extern template class splinepy::splines::Bezier<4, 4>;
extern template class splinepy::splines::Bezier<4, 5>;
extern template class splinepy::splines::Bezier<4, 6>;
extern template class splinepy::splines::Bezier<4, 7>;
extern template class splinepy::splines::Bezier<4, 8>;
extern template class splinepy::splines::Bezier<4, 9>;
extern template class splinepy::splines::Bezier<4, 10>;

extern template class splinepy::splines::Bezier<5, 1>;
extern template class splinepy::splines::Bezier<5, 2>;
extern template class splinepy::splines::Bezier<5, 3>;
extern template class splinepy::splines::Bezier<5, 4>;
extern template class splinepy::splines::Bezier<5, 5>;
extern template class splinepy::splines::Bezier<5, 6>;
extern template class splinepy::splines::Bezier<5, 7>;
extern template class splinepy::splines::Bezier<5, 8>;
extern template class splinepy::splines::Bezier<5, 9>;
extern template class splinepy::splines::Bezier<5, 10>;

extern template class splinepy::splines::Bezier<6, 1>;
extern template class splinepy::splines::Bezier<6, 2>;
extern template class splinepy::splines::Bezier<6, 3>;
extern template class splinepy::splines::Bezier<6, 4>;
extern template class splinepy::splines::Bezier<6, 5>;
extern template class splinepy::splines::Bezier<6, 6>;
extern template class splinepy::splines::Bezier<6, 7>;
extern template class splinepy::splines::Bezier<6, 8>;
extern template class splinepy::splines::Bezier<6, 9>;
extern template class splinepy::splines::Bezier<6, 10>;

extern template class splinepy::splines::Bezier<7, 1>;
extern template class splinepy::splines::Bezier<7, 2>;
extern template class splinepy::splines::Bezier<7, 3>;
extern template class splinepy::splines::Bezier<7, 4>;
extern template class splinepy::splines::Bezier<7, 5>;
extern template class splinepy::splines::Bezier<7, 6>;
extern template class splinepy::splines::Bezier<7, 7>;
extern template class splinepy::splines::Bezier<7, 8>;
extern template class splinepy::splines::Bezier<7, 9>;
extern template class splinepy::splines::Bezier<7, 10>;

extern template class splinepy::splines::Bezier<8, 1>;
extern template class splinepy::splines::Bezier<8, 2>;
extern template class splinepy::splines::Bezier<8, 3>;
extern template class splinepy::splines::Bezier<8, 4>;
extern template class splinepy::splines::Bezier<8, 5>;
extern template class splinepy::splines::Bezier<8, 6>;
extern template class splinepy::splines::Bezier<8, 7>;
extern template class splinepy::splines::Bezier<8, 8>;
extern template class splinepy::splines::Bezier<8, 9>;
extern template class splinepy::splines::Bezier<8, 10>;

extern template class splinepy::splines::Bezier<9, 1>;
extern template class splinepy::splines::Bezier<9, 2>;
extern template class splinepy::splines::Bezier<9, 3>;
extern template class splinepy::splines::Bezier<9, 4>;
extern template class splinepy::splines::Bezier<9, 5>;
extern template class splinepy::splines::Bezier<9, 6>;
extern template class splinepy::splines::Bezier<9, 7>;
extern template class splinepy::splines::Bezier<9, 8>;
extern template class splinepy::splines::Bezier<9, 9>;
extern template class splinepy::splines::Bezier<9, 10>;

extern template class splinepy::splines::Bezier<10, 1>;
extern template class splinepy::splines::Bezier<10, 2>;
extern template class splinepy::splines::Bezier<10, 3>;
extern template class splinepy::splines::Bezier<10, 4>;
extern template class splinepy::splines::Bezier<10, 5>;
extern template class splinepy::splines::Bezier<10, 6>;
extern template class splinepy::splines::Bezier<10, 7>;
extern template class splinepy::splines::Bezier<10, 8>;
extern template class splinepy::splines::Bezier<10, 9>;
extern template class splinepy::splines::Bezier<10, 10>;
#endif
#endif
