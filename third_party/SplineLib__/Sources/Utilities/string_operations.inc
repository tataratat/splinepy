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

template<typename SingleStringOrMultipleStrings>
void Append(String &string, String const &delimiter, SingleStringOrMultipleStrings const string_or_strings) {
  using std::is_same_v;

  if constexpr (std::is_pointer_v<SingleStringOrMultipleStrings> ||  // Character array or
                is_same_v<String, SingleStringOrMultipleStrings>) {  // single string.
    string += (delimiter + string_or_strings);
  } else {
    using ContainedType = typename SingleStringOrMultipleStrings::value_type;
    using std::for_each;

    if constexpr (is_same_v<String, ContainedType>) {  // Standard container of strings.
      for_each(string_or_strings.begin(), string_or_strings.end(),
               [&] (String const &current_string) { string += (delimiter + current_string); });
    } else {  // Nested standard containers of strings.
      for_each(string_or_strings.begin(), string_or_strings.end(),
               [&] (ContainedType const &current_strings) { Append(string, delimiter, current_strings); });
    }
  }
}

template<typename Type>
Type ConvertToNumber(String const &string) {
  using std::is_same_v;

  if constexpr (is_same_v<Type, int>) {
    return std::stoi(string);
  } else if constexpr (is_same_v<Type, long>) {
    return std::stol(string);
  } else if constexpr (is_same_v<Type, long long>) {
    return std::stoll(string);
  } else if constexpr (is_same_v<Type, float>) {
    return std::stof(string);
  } else if constexpr (is_same_v<Type, double>) {
    return std::stod(string);
  } else if constexpr (is_same_v<Type, long double>) {
    return std::stold(string);
  } else if constexpr (is_named_type<Type>) {
    return Type{ConvertToNumber<typename Type::Type_>(string)};
  } else {
    static_assert(error_handling::is_false<Type>,
                  "Only int, long, long long, float, double, long double are supported.");
  }
}

template<typename Type>
Numbers<Type> ConvertToNumbers(String delimited_string, char const &delimiter) {
  StringVector const &strings = SplitAtDelimiter(delimited_string, delimiter);
  Numbers<Type> converted_vector{};
  converted_vector.reserve(strings.size());
  std::for_each(strings.begin(), strings.end(),
                [&] (String const &string) { converted_vector.emplace_back(ConvertToNumber<Type>(string)); });
  return converted_vector;
}

template<typename Type>
String Write(Type const &value, Precision const &precision) {
  static_assert(is_named_type<Type>, "Source type must be a NamedType.");

  std::ostringstream string_stream;
  string_stream.precision(precision.Get());
  string_stream << value.Get();
  return string_stream.str();
}

template<typename ContainerTypeTo, typename ContainerTypeFrom>
ContainerTypeTo Write(ContainerTypeFrom const &from, Precision const &precision) {
  using TypeFrom = typename ContainerTypeFrom::value_type;
  using error_handling::is_false, std::transform, std_container_operations::is_array,
        std_container_operations::is_vector;

  ContainerTypeTo to;
  constexpr bool const is_container_vector = is_vector<ContainerTypeTo>;
  if constexpr (is_array<ContainerTypeTo> || is_container_vector) {
    int const &from_size = from.size();
    if constexpr (is_container_vector) to.reserve(from_size);
    if constexpr (is_named_type<TypeFrom>) {
      std::function<String (TypeFrom const &)> const &write = [&] (TypeFrom const &value) {
                                                                  return Write(value, precision); };
      if constexpr (is_container_vector) {
        transform(from.begin(), from.end(), std::back_inserter(to), write);
      } else {
        transform(from.begin(), from.end(), to.begin(), write);
      }
    } else if constexpr (is_array<TypeFrom> || is_vector<TypeFrom>) {
      using TypeTo = typename ContainerTypeTo::value_type;

      Index::ForEach(0, from_size, [&] (Index const &container) {
          Index::Type_ const &current_container = container.Get();
          if constexpr (is_container_vector) {
            to.emplace_back(Write<TypeTo>(from[current_container], precision));
          } else {
            to[container] = Write<TypeTo>(from[current_container], precision);
          } });
    } else {
      static_assert(is_false<TypeFrom>, "Only NamedTypes and nested std::arrays or std::vectors are supported.");
    }
    return to;
  } else {
    static_assert(is_false<ContainerTypeTo>, "Only std::array and std::vector are supported.");
  }
}
