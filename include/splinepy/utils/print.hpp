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

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace splinepy::utils {

template<typename... Args>
void PrintInfo(Args&&... args) {
  std::cout << "SPLINEPY INFO - ";
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << "\n";
}

/// debug printer - only for debug build
template<typename... Args>
void PrintDebug(Args&&... args) {
#ifndef NDEBUG
  std::cout << "SPLINEPY DEBUG - ";
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << "\n";
#endif
}

template<typename... Args>
void PrintWarning(Args&&... args) {
  std::cout << "SPLINEPY WARNING - ";
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << "\n";
}

template<typename... Args>
void PrintAndThrowError(Args&&... args) {
  std::stringstream error_message{};
  error_message << "SPLINEPY ERROR - ";
  ((error_message << std::forward<Args>(args) << " "), ...);
  error_message << "\n";
  throw std::runtime_error(error_message.str());
}
} /* namespace splinepy::utils */
