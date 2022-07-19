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
    // IRIT expects order not degrees
    export_file << spline.GetDegrees()[i] + 1 << " ";
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
    std::ofstream& export_file, const bool base64encoding) {
  // Increase user readability
  const std::string indentation{"      "};

  // Write Dimensionality
  export_file << indentation
              << "\"dim\" : " << PhysicalPointType::kSpatialDimension << ",\n";

  // Write Parametric Dimensionality
  export_file << indentation << "\"para_dim\" : " << parametric_dimension
              << ",\n";

  // Write degrees
  export_file << indentation << "\"degrees\" : [" << spline.GetDegrees()[0];
  for (std::size_t i{1}; i < parametric_dimension; i++) {
    // IRIT expexts order not degrees
    export_file << ", " << spline.GetDegrees()[i];
  }
  export_file << "],\n";

  // Write All Control Points
  if (base64encoding) {
    export_file << indentation << "\"control_points\" : \""
                << utils::Base64::Encode(spline.control_points) << "\"\n";
  } else {
    export_file << indentation << "\"control_points\" : [\n"
                << indentation << "  " << spline.control_points[0];
    for (std::size_t i_point{1}; i_point < spline.NumberOfControlPoints;
         i_point++) {
      export_file << ",\n"
                  << indentation << "  " << spline.control_points[i_point];
    }
    export_file << "\n" << indentation << "]\n";
  }
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
    export_file << spline.GetDegrees()[i_par_dim]
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
    for (std::size_t i_point{0};
         i_point < 2 * spline.GetDegrees()[i_par_dim] + 2; i_point++) {
      export_file << indentation
                  << (i_point < (spline.GetDegrees()[i_par_dim] + 1) ? "0."
                                                                     : "1.");
    }
    export_file << "\n" << indentation << indentation << "</kntVec>\n";
  }
  export_file << indentation << "</kntVecs>\n";

  // Write Footer
  export_file << indentation << "</SplineEntry>\n";
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::GuessByExtension(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  GuessByExtension(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::GuessByExtension(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // In case of a absolute path strip filename
  std::size_t last_found_id = filename.find_last_of("/\\");
  const std::string file_without_path = filename.substr(last_found_id + 1);

  // Check file extension
  std::size_t extension_start_id = file_without_path.find_last_of(".");
  const std::string extension_str =
      file_without_path.substr(extension_start_id + 1);

  if (extension_str == "itd") {
    Logger::Logging("Identified IRIT export format");
    // Export as IRIT
    AsIRIT(spline_group, filename);
  } else if (extension_str == "xml") {
    Logger::Logging("Identified XML export format");
    // Export as XML
    AsXML(spline_group, filename);
  } else if (extension_str == "json") {
    Logger::Logging("Identified json export format");
    Logger::ExtendedInformation("Json Export format uses base64 as default");
    Logger::ExtendedInformation("Use AsJSON function to change to text output");
    // Export as custom JSON
    AsJSON(spline_group, filename);
  } else if (extension_str == "mesh") {
    Logger::Logging("Identified MFEM export format");
    // Export as MFEM
    if constexpr (parametric_dimension ==
                      PhysicalPointType::kSpatialDimension &&
                  (parametric_dimension == static_cast<std::size_t>(2) ||
                   parametric_dimension == static_cast<std::size_t>(3))) {
      AsMFEM(spline_group, filename);
    } else {
      Logger::Error("MFEM only supports 2D->2D or 3D->3D Splines");
      Logger::TerminatingError("Not Implemented");
    }
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
  Logger::UserInfo("Starting IRIT Export as : " + filename);
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
    // Export individual splines separately
    format2IRITfile(spline_group[i], export_file);
    // Finalize Current Spline
    export_file << "\n\t]";
  }

  // Write Footer
  export_file << "\n]";
  export_file.close();
  Logger::UserInfo("Finished IRIT Export");
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
  Logger::UserInfo("Starting XML Export as : " + filename);
  Logger::Warning(
      "XML Export inefficient (export as NURBS) and not standardized.");
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
    // Export individual splines separately
    format2XMLfile(spline_group[i], export_file);
  }

  // Write Footer
  export_file << "</SplineList>\n";
  export_file.close();
  Logger::UserInfo("Finished XML Export");
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename, const bool base64encoding) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsJSON(spline_group, filename, base64encoding);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename, const bool base64encoding) {
  Logger::UserInfo("Starting JSON Export as : " + filename);
  Logger::Logging(base64encoding ? "Exporting using Base64 Encoding"
                                 : "Exporting as Text-Only");
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "json"));

  // Write Header
  export_file << "{\n  \"Name\" : \"BezierSplineGroup\",\n";

  // Number of Splines
  export_file << "  \"NumberOfSplines\" : " << spline_group.size() << ",\n";

  // Encoding
  export_file << "  \"Base64Encoding\" : " << base64encoding << ",\n";

  // Start Spline Export
  export_file << "  \"SplineList\" : [\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // ID and Type are set here
    export_file
        << "    {\n      \"SplineType\" : \"Bezier\",\n      \"SplineID\" : "
        << i << ",\n";
    // Export individual splines separately
    format2JSONfile(spline_group[i], export_file, base64encoding);
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
  Logger::UserInfo("Finished XML Export");
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsMFEM(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsMFEM(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsMFEM(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  Logger::UserInfo("Starting MFEM Export as : " + filename);
  Logger::Warning(
      "MFEM Exports of derivatives will result in different enumerations!");
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "mesh"));

  // MFEM format currently restricted to 2D2D and 3D3D
  static_assert(parametric_dimension == PhysicalPointType::kSpatialDimension,
                "Mismatch between parametric and physical dimension");
  static_assert(parametric_dimension == static_cast<std::size_t>(2) ||
                    parametric_dimension == static_cast<std::size_t>(3),
                "MFEM export limited to 2D and 3D meshes");

  // -- Auxiliary values --
  const std::size_t number_of_splines{spline_group.size()};
  constexpr auto subelement_vertex_ids =
      algorithms::HyperCube<parametric_dimension>::SubElementVerticesToFace();
  constexpr std::size_t number_of_relevant_element_faces{
      2 * parametric_dimension - 2};
  constexpr std::size_t number_of_vertices_per_subelement{algorithms::IntPower(
      static_cast<std::size_t>(2), parametric_dimension - 1)};
  constexpr std::size_t number_of_element_vertices{
      algorithms::IntPower(static_cast<std::size_t>(2), parametric_dimension)};
  constexpr auto faces_orthogonal_to_parametric_dimension =
      algorithms::HyperCube<
          parametric_dimension>::GetNormalFaceIndicesToParametricDimension();

  // Retrieve Connectivity
  const auto connectivity =
      bezman::utils::algorithms::GetConnectivityForSplineGroup(
          spline_group);
  if (connectivity.size() == 0) {
    Logger::TerminatingError("Creating connectivity failed.");
  }
  Logger::ExtendedInformation("Connectivity retrieved without errors");

  // Initialize counters
  std::size_t number_of_boundaries{};

  // Initialize Vertices to identify their vertex ids
  std::vector<PhysicalPointType> corner_vertices;
  corner_vertices.reserve(number_of_splines * number_of_element_vertices);

  // Retrieve all vertices
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    const auto global_vertex_id =
        algorithms::HyperCube<parametric_dimension>::VertexIdForDegrees(
            spline_group[i_spline].GetDegrees());
    for (std::size_t i_point{}; i_point < number_of_element_vertices;
         i_point++) {
      corner_vertices.push_back(
          spline_group[i_spline].control_points[global_vertex_id[i_point]]);
    }
  }
  // Retrieve Vertex enumeration
  const auto vertex_id = algorithms::IndexUniquePointList(
      corner_vertices,
      spline_group.MaximumCorner() - spline_group.MinimumCorner());
  Logger::ExtendedInformation("Enumerated Vertices using uniquified points");

  // Continue with 2.
  // Enumerate knot vectors
  std::vector<std::array<std::size_t, parametric_dimension>> knot_vector_ids(
      number_of_splines,
      // Fill with array created by lambda function
      []() {
        std::array<std::size_t, parametric_dimension> init_{};
        init_.fill(static_cast<std::size_t>(-1));
        return init_;
      }());

  // Initialize knot_vector_counter
  std::size_t n_assigned_knot_vectors{};

  // index of splines which has to be control if they have neighbors
  std::vector<std::size_t> splines_to_go{};

  // loop over all splines
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    // loop over all dimensions
    for (std::size_t i_dim{}; i_dim < parametric_dimension; i_dim++) {
      // check if spline has already a knot vector id
      if (knot_vector_ids[i_spline][i_dim] != static_cast<size_t>(-1)) {
        continue;
      }

      // insert index of start spline
      splines_to_go.push_back(i_spline);

      // retrieve relevant faces for propagation
      const auto& faces = faces_orthogonal_to_parametric_dimension[i_dim];

      // Assign knot vector id to current spline
      knot_vector_ids[i_spline][i_dim] = n_assigned_knot_vectors;

      // run through all splines until the vector splines_to_go is empty
      while (!splines_to_go.empty()) {
        // get last spline_id of the vector
        const std::size_t current_spline_id = splines_to_go.back();
        // delete the set current_spline_id from the vector
        splines_to_go.pop_back();

        // run through all neighbors
        for (std::size_t i_neighbor{};
             i_neighbor < number_of_relevant_element_faces; i_neighbor++) {
          const std::size_t& face_to_check = faces[i_neighbor];
          // Check if boundary
          if (connectivity[current_spline_id][face_to_check] ==
              static_cast<std::size_t>(-1)) {
            number_of_boundaries++;
            continue;
          }
          const std::size_t& current_neighbor =
              connectivity[current_spline_id][face_to_check];

          // if neighbor spline has the same  knot vector id as the
          // current spline continue --> Loop is closed
          if (knot_vector_ids[current_neighbor][i_dim] ==
              n_assigned_knot_vectors) {
            continue;
          }

          // Check if neighbor is different knot_vector_id, if true something
          // went wrong --> error
          if (knot_vector_ids[current_neighbor][i_dim] !=
              static_cast<std::size_t>(-1)) {
            Logger::TerminatingError(
                "Connectivity error. Knotvector mismatch between two "
                "neighboring splines");
          }
          if (spline_group[current_neighbor].GetDegrees()[i_dim] !=
              spline_group[i_spline].GetDegrees()[i_dim]) {
            Logger::TerminatingError(
                "Degree mismatch between two neighboring splines");
          }

          // set new knot vector
          knot_vector_ids[current_neighbor][i_dim] = n_assigned_knot_vectors;
          // Add neighbor spline to splines_to_go
          splines_to_go.push_back(current_neighbor);
        }
      }
      n_assigned_knot_vectors++;
    }
  }
  Logger::ExtendedInformation("Knotvector assignment successful");

  // The boundary count is only valid for 2D elements, as all faces correspond
  // to exactly one knot_vector dimension. In higher order problems, like 3D the
  // same face needs to be checked for n_dim-1 knot_vector dimensions, the
  // number of boundaries therefore needs to be updated accordingly
  assert(number_of_boundaries % (parametric_dimension - 1) == 0);
  number_of_boundaries /= (parametric_dimension - 1);

  // map for the edges -> key is the edge with the vertices | value is the
  // knot vector id
  std::map<std::pair<std::size_t, std::size_t>, std::size_t>
      edge_uniquify_map{};

  constexpr auto local_face_vertex_index =
      algorithms::HyperCube<parametric_dimension>::EdgeVertexIndices();
  constexpr std::size_t number_of_edges_per_element =
      local_face_vertex_index.size();

  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    for (std::size_t i_edge{}; i_edge < number_of_edges_per_element; i_edge++) {
      // Edges are ordered by dimension using integer division
      std::size_t i_dim{i_edge / (parametric_dimension == 2 ? 2 : 4)};
      // Add element to edge
      edge_uniquify_map[
          // Create Pair of vertex ids as key
          std::make_pair(
              // First Vertex
              vertex_id[i_spline * number_of_element_vertices +
                        local_face_vertex_index[i_edge][0]],
              // Second Vertex
              vertex_id[i_spline * number_of_element_vertices +
                        local_face_vertex_index[i_edge][1]])] =
          // set the value
          knot_vector_ids[i_spline][i_dim];
    }
  }

  // Write Header
  export_file << "MFEM NURBS mesh v1.0\n";

  // Write dimension
  export_file << "\ndimension \n" << parametric_dimension << "\n";

  // Write segment element
  export_file << "\nelements\n";
  export_file << number_of_splines << "\n";

  // Maximum vertex_id for vertex_output
  std::size_t maximum_vertex_id{0};

  // Write Elements
  // structure of line:  <element attribute><geometric type><vertex 0>...
  // <vertex m>
  // run through all splines
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    // write the first two entries of the line:
    // geometric type: 3 = square | 5 = cube
    export_file << 1 << " " << ((parametric_dimension == 2) ? 3 : 5);
    // run through all vertex inside the current line
    for (std::size_t i_vertex{}; i_vertex < number_of_element_vertices;
         i_vertex++) {
      const std::size_t& c_vertex =
          vertex_id[i_spline * number_of_element_vertices + i_vertex];
      export_file << " " << c_vertex;
      // Update maximum vertex_id, this is needed inside the vertices segment
      maximum_vertex_id =
          (c_vertex > maximum_vertex_id) ? c_vertex : maximum_vertex_id;
    }
    export_file << "\n";
  }

  // Write out boundaries
  // structure of line:  <boundary id><boundary element type><vertex 0>...
  // <vertex m>
  export_file << "\nboundary\n";
  // write the number of boundaries
  export_file << number_of_boundaries << "\n";
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    for (std::size_t i_face{}; i_face < 2 * parametric_dimension; i_face++) {
      if (connectivity[i_spline][i_face] == static_cast<std::size_t>(-1)) {
        // write the first two entries of the line:
        // geometric type: 1 = line | 3 = square
        export_file << 1 << " " << ((parametric_dimension == 2) ? 1 : 3);
        // write the vertices
        const auto& local_face_vertices = subelement_vertex_ids[i_face];
        for (std::size_t i_local_vertex_on_face{};
             i_local_vertex_on_face < local_face_vertices.size();
             i_local_vertex_on_face++) {
          export_file << " "
                      << vertex_id[i_spline * number_of_element_vertices +
                                   local_face_vertices[i_local_vertex_on_face]];
        }
        export_file << "\n";
      }
    }
  }

  // Export Edges (Corner line segments)
  // structure of line:  <knot vector id><vertex start><vertex end>
  export_file << "\n\nedges\n";
  // get the number of edges from the map
  const std::size_t n_edges = edge_uniquify_map.size();
  // write the number of edges
  export_file << n_edges << "\n";
  for (const auto& [edge_vertices, i_knot_id] : edge_uniquify_map) {
    // write edges
    export_file << i_knot_id << " " << edge_vertices.first << " "
                << edge_vertices.second << "\n";
  }
  export_file << "\n";

  // Export Number of Vertices
  export_file << "\n\nvertices\n";
  // write the number of vertices
  export_file << maximum_vertex_id + 1 << "\n";

  using BezierSplineType =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;

  // Export all Splines as patches and their control points
  export_file << "\n\npatches\n";
  // run through all splines/patches and write the knot vector- and control
  // point segments
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    export_file << "\nknotvectors"
                << "\n";
    // write the number of knot vectors
    export_file << BezierSplineType::kParametricDimensions << "\n";

    // structure of line <degree><number of control points><xi 0>...<xi m>
    // copy code from format2XMLfile  ;)
    for (std::size_t i_par_dim{};
         i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
      // Write All Control Points
      export_file << spline_group[i_spline].GetDegrees()[i_par_dim]
                  << " "
                  // Number of control points in parametric dimension
                  << spline_group[i_spline].GetDegrees()[i_par_dim] + 1 << " ";
      for (std::size_t i_point{};
           i_point < 2 * spline_group[i_spline].GetDegrees()[i_par_dim] + 2;
           i_point++) {
        export_file
            << (i_point < (spline_group[i_spline].GetDegrees()[i_par_dim] + 1)
                    ? "0.0 "
                    : "1.0 ");
      }
      export_file << "\n";
    }

    // write physical dimension
    export_file << "\ndimension"
                << "\n";
    export_file << BezierSplineType::kParametricDimensions << "\n";

    // write the control point segment
    // structure of line <x><y>(<z>)<weight>
    export_file << "\n\ncontrolpoints_cartesian\n";
    for (std::size_t i_point{};
         i_point < spline_group[i_spline].NumberOfControlPoints; i_point++) {
      for (std::size_t i_dim{}; i_dim < PhysicalPointType ::kSpatialDimension;
           i_dim++) {
        export_file << spline_group[i_spline].control_points[i_point][i_dim]
                    << " ";
      }
      // add weights (polynomial splines are always unweighted)
      export_file << 1. << "\n";
    }

    export_file << "\n";
  }

  // Write Footer
  export_file.close();
  Logger::UserInfo("Finised XML Export");
}

#endif  // SRC_UTILS_EXPORT_INC
