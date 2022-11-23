#pragma once

#include <iostream>
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
  std::cout << "SPLINEPY ERROR - ";
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << "\n";
  throw std::runtime_error("Error Occured! Abort the mission!");
}
} /* namespace splinepy::utils */
