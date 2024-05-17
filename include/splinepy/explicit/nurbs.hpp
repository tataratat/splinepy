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
extern template class splinepy::splines::Nurbs<1>;
extern template class splinepy::splines::Nurbs<2>;
extern template class splinepy::splines::Nurbs<3>;
#ifdef SPLINEPY_MORE
extern template class splinepy::splines::Nurbs<4>;
extern template class splinepy::splines::Nurbs<5>;
extern template class splinepy::splines::Nurbs<6>;
extern template class splinepy::splines::Nurbs<7>;
extern template class splinepy::splines::Nurbs<8>;
extern template class splinepy::splines::Nurbs<9>;
extern template class splinepy::splines::Nurbs<10>;
#endif
#endif
