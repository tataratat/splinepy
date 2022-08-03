/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

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

#ifndef SRC_UTILS_LOGGER_HPP
#define SRC_UTILS_LOGGER_HPP

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace bezman::utils {
/**
 * @brief Provides functionality for output information
 *
 * The logger class provides functionality to write log files and to format user
 * information, warnings and errors for the console output. There only ever
 * exists a single instance of the Logger class, as it is provided as a
 * singleton object. The Level of information that is provided by the logger can
 * be changed at runtime, however, the logger needs to be activated with a
 * Preprocessor flag at compile time
 */
class Logger {
 public:
  /**
   * @brief Describes the scope of output information
   *
   * Each enum object describes a specific bit of information and if it is to be
   * respected or not
   */
  enum class OutputLevel : unsigned int {
    nothing = 0,
    userinfo = 1,
    errors = 2,
    warnings = 4,
    logging = 8,
    logging_verbose = 16,
    time_stamp = 32,
    warning_file = 64,
    error_file = 128,
    log_file = 256,
    all = 511
  };

  /**
   * @brief Set the Output Level
   *
   * Set the output level, default its just written into the terminal. The
   * outputlevel is set via an unsigned integer value which represents the sum
   * of all chosen enum objects, see OutputLevel
   */
  static void SetOutputLevel(unsigned int outputlevel) {
    Get().SetOutputLevel_(outputlevel);
  }

  /**
   * @brief Set the Output Level with initializer list
   *
   * Set the output level, default its just written into the terminal. The
   * outputlevel is set via a list which contains all chosen enum objects,
   * see OutputLevel
   */
  static void SetOutputLevel(
      const std::initializer_list<OutputLevel>& output_options) {
    Get().SetOutputLevel_(output_options);
  }

  /**
   * @brief Write formatted error warnings into the terminal and the error file
   *
   * This function writes formatted output, but also throws an exception which
   * potentially leads to the termination of the program
   */
  static void Warning(const std::string& warning_text) {
    Get().Warning_(warning_text);
  }

  /**
   * @brief Write formatted error warnings into the terminal and the error file
   *
   * This function writes formatted output, but also throws an exception which
   * potentially leads to the termination of the program
   */
  template <typename ExceptionType = std::runtime_error>
  static void TerminatingError(const std::string& error_text) {
    Get().TerminatingError_<ExceptionType>(error_text);
  }

  /**
   * @brief Write formatted error warnings into the terminal and the error file
   *
   * This function writes formatted output, but also terminates the execution of
   * the problem
   */
  static void Error(const std::string& error_text) { Get().Error_(error_text); }

  /**
   * @brief Write a formatted user information into the terminal and log file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  static void UserInfo(const std::string& info_text) {
    Get().UserInfo_(info_text);
  }

  /**
   * @brief Write a formatted logging information into the terminal and log file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  static void Logging(const std::string& log_text) { Get().Logging_(log_text); }

  /**
   * @brief Write a formatted extended information into the terminal and log
   * file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  static void ExtendedInformation(const std::string& log_text) {
    Get().ExtendedInformation_(log_text);
  }

 private:
  /**
   * @brief Get the Logger object
   *
   * Returns a singleton reference to the logger instance.
   */
  static Logger& Get() {
    static Logger singleton_instance;
    singleton_instance.init();
    return singleton_instance;
  }

  /**
   * @brief Set the Output Level
   *
   * Set the output level, default its just written into the terminal. The
   * outputlevel is set via an unsigned integer value which represents the sum
   * of all chosen enum objects, see OutputLevel
   */
  void SetOutputLevel_(unsigned int outputlevel) {
    outputlevel_ = outputlevel;
    init();
  }

  /**
   * @brief Set the Output Level
   *
   * Set the output level, default its just written into the terminal. The
   * outputlevel is set via a series of OutputLevel options, possible call is
   *
   * SetOutputLevel({OutputLevel::warnings, OutputLevel::errors});
   */
  void SetOutputLevel_(
      const std::initializer_list<OutputLevel>& output_options) {
    outputlevel_ = 0;
    for (auto i_option = output_options.begin();
         i_option != output_options.end(); ++i_option) {
      outputlevel_ |= static_cast<unsigned>(*i_option);
    }
    init();
  }

  /**
   * @brief Write a formatted warning into the terminal and warning log file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  void Warning_(const std::string& warning_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::warnings) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("warning", Color::yellow) << "]" << GetTimeStamp_()
                << " : " << ColorText(warning_text, Color::yellow) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      warning_file << "[" << std::setw(padding_first_col_files) << "warning"
                   << "]" << GetTimeStamp_() << " : " << warning_text << "\n";
    }
#endif
  }

  /**
   * @brief Write a formatted user information into the terminal and log file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  void UserInfo_(const std::string& info_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::userinfo) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("user info", Color::blue) << "]" << GetTimeStamp_()
                << " : " << ColorText(info_text, Color::blue) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_col_files) << "user info"
               << "]" << GetTimeStamp_() << " : " << info_text << "\n";
    }
#endif
  }

  /**
   * @brief Write a formatted logging information into the terminal and log file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  void Logging_(const std::string& log_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::logging) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("log", Color::shell_default) << "]"
                << GetTimeStamp_() << " : "
                << ColorText(log_text, Color::shell_default) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_col_files) << "user info"
               << "]" << GetTimeStamp_() << " : " << log_text << "\n";
    }
#endif
  }

  /**
   * @brief Write a formatted extended information into the terminal and log
   * file
   *
   * The function writes formatted output based on the chosen level of
   * information
   */
  void ExtendedInformation_(const std::string& log_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::logging_verbose) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("log", Color::green) << "]" << GetTimeStamp_()
                << " : " << ColorText(log_text, Color::green) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_col_files) << "user info"
               << "]" << GetTimeStamp_() << " : " << log_text << "\n";
    }
#endif
  }

  /**
   * @brief Write formatted error warnings into the terminal and the error file
   *
   * This function writes formatted output, but also terminates the execution of
   * the problem
   */
  void Error_(const std::string& error_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::errors) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("error", Color::red) << "]" << GetTimeStamp_()
                << " : " << ColorText(error_text, Color::red) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      error_file << "[" << std::setw(padding_first_col_files) << "error"
                 << "]" << GetTimeStamp_() << " : " << error_text << "\n";
    }
#endif
  }

  /**
   * @brief Write formatted error warnings into the terminal and the error file
   *
   * This function writes formatted output, but also throws an exception which
   * potentially leads to the termination of the program
   */
  template <typename ExceptionType = std::runtime_error>
  void TerminatingError_(const std::string& error_text) {
    if (static_cast<unsigned>(OutputLevel::errors) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_col_console)
                << ColorText("crit. error", Color::red) << "]"
                << GetTimeStamp_() << " : " << ColorText(error_text, Color::red)
                << "\n";
    }
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      error_file << "[" << std::setw(padding_first_col_files) << "error"
                 << "]" << GetTimeStamp_() << " : " << error_text << "\n";
    }
#endif
    // Close file streams in order to save error files
    close();
    // Throw exception
    throw ExceptionType(error_text);
  }

  /**
   * @brief Destroy the Logger object and close files
   */
  ~Logger() { close(); }

 private:
  /// Padding in console
  const unsigned int padding_first_col_console{20};

  /// Padding in text files
  const unsigned int padding_first_col_files{10};

  /// File Streams
  std::ofstream log_file;
  std::ofstream warning_file;
  std::ofstream error_file;

  /// Initialize output streams
  void init() {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      if (!log_file.is_open()) {
        log_file.open("log_file.log");
      }
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      if (!warning_file.is_open()) {
        warning_file.open("warning_file.log");
      }
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      if (!error_file.is_open()) {
        error_file.open("error_file.log");
      }
    }
#endif
  }

  /// Close File-Streams to save log files
  void close() {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      if (log_file.is_open()) {
        log_file.close();
      }
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      if (warning_file.is_open()) {
        warning_file.close();
      }
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      if (error_file.is_open()) {
        error_file.close();
      }
    }
#endif
  }

  /// Output Level 39 {: user_info, error, warning with timestamps}
  unsigned int outputlevel_{39};

  /**
   * @brief Get the Time Stamp as a string in brackets
   */
  std::string GetTimeStamp_() {
    if (static_cast<unsigned>(OutputLevel::time_stamp) & outputlevel_) {
      std::time_t result = std::time(nullptr);
      std::string time_stamp = (std::asctime(std::localtime(&result)));
      time_stamp.pop_back();
      return "[" + time_stamp + "]";
    } else {
      return std::string();
    }
  }

  // Available Colors
  enum Color { red, green, yellow, blue, white, cyan, magenta, shell_default };

  /// Color Text
  std::string ColorText(const std::string& text, Color color) {
    switch (color) {
      case Color::red:
        return std::string("\033[31m") + text + std::string("\033[0m");
      case Color::green:
        return std::string("\033[32m") + text + std::string("\033[0m");
      case Color::yellow:
        return std::string("\033[33m") + text + std::string("\033[0m");
      case Color::blue:
        return std::string("\033[34m") + text + std::string("\033[0m");
      case Color::magenta:
        return std::string("\033[35m") + text + std::string("\033[0m");
      case Color::cyan:
        return std::string("\033[36m") + text + std::string("\033[0m");
      case Color::white:
        return std::string("\033[37m") + text + std::string("\033[0m");
      case Color::shell_default:
        return std::string("\033[;0m") + text + std::string("\033[0m");
        return text;
      default:
        return text;
    }
  }

  // Prohibit implementation of Loggers constructors by making them private
  Logger() {}
  Logger(Logger const&);
  void operator=(Logger const&);
};
}  // namespace bezman::utils
#endif  // SRC_UTILS_LOGGER_HPP
