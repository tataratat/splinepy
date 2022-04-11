#ifndef SRC_UTILS_EXPORT_INC
#define SRC_UTILS_EXPORT_INC

std::string Export::ensureFileExtension(const std::string& filename,
                                        const std::string& obj_extension) {
  std::size_t ext_start = filename.find_last_of(".");
  const std::string extension = filename.substr(ext_start + 1);
  if (extension == obj_extension) {
    return filename;
  } else {
    return (filename + ".") + obj_extension;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2IRITfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file) {
  // Irit only support 9 dimensional spaces
  static_assert((PhysicalPointType::kSpatialDimension < 10),
                "IRIT natively only supports up to 9D spaces");

  // Increase user readability
  const std::string indentation{"\t"};

  // Write Header
  export_file << indentation << indentation << "[MULTIVAR BEZIER "
              << parametric_dimension << " ";

  // Write different degrees
  for (std::size_t i{}; i < parametric_dimension; i++) {
    // IRIT expexts order not degrees
    export_file << spline.degrees[i] + 1 << " ";
  }

  // Point Type
  export_file << "E" << PhysicalPointType::kSpatialDimension << "\n";

  // Write All Control Points
  for (std::size_t i_point{0}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    export_file << indentation << indentation << indentation
                << spline.control_points[i_point] << "\n";
  }

  // Write Footer
  export_file << indentation << indentation << "]";
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2JSONfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file) {
  // Increase user readability
  const std::string indentation{"      "};

  // Write Dimensionality
  export_file << indentation
              << "\"Dimension\" : " << PhysicalPointType::kSpatialDimension
              << ",\n";

  // Write Parametric Dimensionality
  export_file << indentation
              << "\"ParametricDimension\" : " << parametric_dimension << ",\n";

  // Write degrees
  export_file << indentation << "\"Degrees\" : [" << spline.degrees[0];
  for (std::size_t i{1}; i < parametric_dimension; i++) {
    // IRIT expexts order not degrees
    export_file << ", " << spline.degrees[i];
  }
  export_file << "],\n";

  // Write All Control Points
  export_file << indentation << "\"ControlPoints\" : [\n"
              << indentation << "  " << spline.control_points[0];
  for (std::size_t i_point{0}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    export_file << ",\n"
                << indentation << "  " << spline.control_points[i_point];
  }
  export_file << "\n" << indentation << "]\n";
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2XMLfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file) {
  // Handy alias
  using BezierSplineType =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;
  // Increase user readability
  const std::string indentation{"  "};

  // Write Header
  export_file << indentation << "<SplineEntry splDim=\"" << parametric_dimension
              << "\" spaceDim=\"" << PhysicalPointType::kSpatialDimension
              << "\" numOfCntrlPntVars=\""
              << PhysicalPointType ::kSpatialDimension << "\" numCntrlPnts=\""
              << spline.NumberOfControlPoints << "\" numOfEleVars=\"" << 0
              << "\">\n";

  // Write ControlPoint Variable Names
  export_file << indentation << "<cntrlPntVarNames>\n";
  for (std::size_t i{}; i < PhysicalPointType ::kSpatialDimension; i++) {
    export_file << indentation << "  x" << i << " ";
  }
  export_file << "\n" << indentation << "</cntrlPntVarNames>\n";

  // Write Points
  export_file << indentation << "<cntrlPntVars>";
  for (std::size_t i_point{}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    // Spacing
    export_file << "\n" << indentation << indentation;
    for (std::size_t i_dim{}; i_dim < PhysicalPointType ::kSpatialDimension;
         i_dim++) {
      export_file << std::setw(12) << std::setprecision(10)
                  << spline.control_points[i_point][i_dim] << indentation;
    }
  }
  export_file << "\n" << indentation << "</cntrlPntVars>\n";

  // Write Weights
  export_file << indentation << "<wght>\n" << indentation << indentation;
  for (std::size_t i_point{}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    export_file << "1.0  ";  // Add \n if required
  }
  export_file << "\n" << indentation << "</wght>\n";

  // Write Degrees
  export_file << indentation << "<deg>\n" << indentation << indentation;
  for (std::size_t i_par_dim{};
       i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
    export_file << spline.degrees[i_par_dim]
                << indentation;  // Add \n if required
  }
  export_file << "\n" << indentation << "</deg>\n";

  // Knotvectors (required)
  export_file << indentation << "<kntVecs>\n";
  for (std::size_t i_par_dim{};
       i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
    export_file << indentation << indentation << "<kntVec>\n"
                << indentation << indentation;
    // Write All Control Points
    for (std::size_t i_point{0}; i_point < 2 * spline.degrees[i_par_dim] + 2;
         i_point++) {
      export_file << indentation
                  << (i_point < (spline.degrees[i_par_dim] + 1) ? "0." : "1.");
    }
    export_file << "\n" << indentation << indentation << "</kntVec>\n";
  }
  export_file << indentation << "</kntVecs>\n";

  // Write Footer
  export_file << indentation << "</SplineEntry>\n";
}

template <typename ExportType>
void Export::GuessByExtension(const ExportType& export_object,
                              const std::string& filename) {
  // In case of a absolute path strip filename
  std::size_t last_found_id = filename.find_last_of("/\\");
  const std::string file_without_path = filename.substr(last_found_id + 1);

  // Check file extension
  std::size_t extension_start_id = file_without_path.find_last_of(".");
  const std::string extension_str =
      file_without_path.substr(extension_start_id + 1);

  if (extension_str == "itd") {
    // Export as IRIT
    AsIRIT(export_object, filename);
  } else if (extension_str == "xml") {
    // Export as XML
    AsXML(export_object, filename);
  } else if (extension_str == "json") {
    // Export as custom JSON
    AsJSON(export_object, filename);
  } else {
    throw std::invalid_argument("File Extension not supported");
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsIRIT(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsIRIT(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsIRIT(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "itd"));

  // Write Header
  export_file << "[OBJECT BEZIERGROUP";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // Write Header
    export_file << "\n\t[OBJECT BEZIERSPLINE" << i;
    // Small spacer for Readability
    export_file << "\n";
    // Export individual splines seperately
    format2IRITfile(spline_group[i], export_file);
    // Finalize Current Spline
    export_file << "\n\t]";
  }

  // Write Footer
  export_file << "\n]";
  export_file.close();
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsXML(const BezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& spline,
                   const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsXML(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsXML(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "xml"));

  // Write Header
  export_file << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";

  // Start Spline List
  export_file << "<SplineList SplineType = \"1\" NumberOfSplines =\""
              << spline_group.size() << "\">\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // Export individual splines seperately
    format2XMLfile(spline_group[i], export_file);
  }

  // Write Footer
  export_file << "</SplineList>\n";
  export_file.close();
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsXML(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "json"));

  // Write Header
  export_file << "{\n  \"Name\" : \"BezierSplineGroup\",\n";

  // Number of Splines
  export_file << "  \"NumberOfSplines\" : " << spline_group.size() << ",\n";

  // Start Spline Export
  export_file << "  \"SplineList\" : [\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // ID and Type are set here
    export_file
        << "    {\n      \"SplineType\" : \"Bezier\",\n      \"SplineID\" : "
        << i << ",\n";
    // Export individual splines seperately
    format2JSONfile(spline_group[i], export_file);
    export_file << "    }";
    if (i < spline_group.size() - 1) {
      export_file << ",\n";
    } else {
      export_file << "\n";
    }
  }

  // Write Footer
  export_file << "  ]\n}";
  export_file.close();
}

#endif  // SRC_UTILS_EXPORT_INC
