/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

#include "Sources/Utilities/string_operations.hpp"

#include <algorithm>
#include <functional>
#include <sstream>
#include <utility>

namespace splinelib::sources::utilities::string_operations {

bool EndsWith(String const &string, String const &pattern) {
  int const &length_of_string = string.length(), &length_of_pattern = pattern.length();
  return length_of_string >= length_of_pattern ? string.compare((length_of_string - length_of_pattern),
                                                                length_of_pattern, pattern) == 0 : false;
}

bool StartsWith(String const &string, String const &pattern) {
  size_t const &length_of_pattern = pattern.length();
  return string.length() >= length_of_pattern ? string.compare(0, length_of_pattern, pattern) == 0 : false;
}

String TrimCharacter(String string, char const &character) {
  using std::find_if;

  String::iterator const &string_begin = string.begin();
  std::function<bool(char const &)> const &is_not_character = [&] (char const &current_character) {
                                                                  return (current_character != character); };
  string.erase(string_begin, find_if(string_begin, string.end(), is_not_character));
  string.erase(find_if(string.rbegin(), string.rend(), is_not_character).base(), string.end());
  return string;
}

StringVector SplitAtDelimiter(String string, char const &delimiter) {
  using StringStream = std::stringstream;
  using std::getline;

  StringStream stringstream_from_string{std::move(string)};
  String current_line;
  StringVector converted_vector;
  while (getline(stringstream_from_string, current_line)) {
    StringStream stringstream_from_line{current_line};
    String current_string;
    while (getline(stringstream_from_line, current_string, delimiter))
        if (!current_string.empty()) converted_vector.emplace_back(current_string);
  }
  return converted_vector;
}

}  // namespace splinelib::sources::utilities::string_operations
