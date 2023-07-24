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

/// debug printer - first argument is bool, so <on, off> is switchable.
template<typename... Args>
void PrintDebug(bool on, Args&&... args) {
  if (on) {
    std::cout << "SPLINEPY DEBUG - ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << "\n";
  }
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
