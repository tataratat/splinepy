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

#include "Sources/InputOutput/operations.hpp"

#include "Sources/InputOutput/iges.hpp"
#include "Sources/InputOutput/irit.hpp"
#include "Sources/InputOutput/xml.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/string_operations.hpp"

namespace splinelib::sources::input_output::operations {

FileFormat DetermineFileFormat(String const &file_name) {
  using utilities::string_operations::EndsWith;

  if (EndsWith(file_name, ".iges")) {
    return FileFormat::kIges;
  } else if (EndsWith(file_name, ".itd")) {
    return FileFormat::kIrit;
  } else if (EndsWith(file_name, ".vtk")) {
    return FileFormat::kVtk;
  } else if (EndsWith(file_name, ".xml")) {
    return FileFormat::kXml;
  } else {
#ifndef NDEBUG
    throw RuntimeError("Only IGES (.iges), IRIT (.itd), VTK (.vtk), and XML (.xml) file formats are supported.");
#endif
    return FileFormat::kInvalid;
  }
}

Splines Read(String const &file_name) {
  switch (DetermineFileFormat(file_name)) {
    case FileFormat::kIges:
      return iges::Read(file_name);
      break;
    case FileFormat::kIrit:
      return irit::Read(file_name);
      break;
    case FileFormat::kVtk:
#ifndef NDEBUG
      throw RuntimeError("VTK (.vtk) files cannot be read.");
#endif
      return Splines{};
      break;
    case FileFormat::kXml:
      return xml::Read(file_name);
      break;
    default:  // case FileFormat::kInvalid
      return Splines{};
      break;
  }
}

}  // namespace splinelib::sources::input_output::operations
