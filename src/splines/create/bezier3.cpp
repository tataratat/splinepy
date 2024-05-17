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

#include <splinepy/splines/create/create_bezier.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier3(const int dim, const int* degrees, const double* control_points) {
  switch (dim) {
  case 1:
    return std::make_shared<Bezier<3, 1>>(degrees, control_points);
  case 2:
    return std::make_shared<Bezier<3, 2>>(degrees, control_points);
  case 3:
    return std::make_shared<Bezier<3, 3>>(degrees, control_points);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<Bezier<3, 4>>(degrees, control_points);
  case 5:
    return std::make_shared<Bezier<3, 5>>(degrees, control_points);
  case 6:
    return std::make_shared<Bezier<3, 6>>(degrees, control_points);
  case 7:
    return std::make_shared<Bezier<3, 7>>(degrees, control_points);
  case 8:
    return std::make_shared<Bezier<3, 8>>(degrees, control_points);
  case 9:
    return std::make_shared<Bezier<3, 9>>(degrees, control_points);
  case 10:
    return std::make_shared<Bezier<3, 10>>(degrees, control_points);
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBezier. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    break;
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateBezier. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines::create
