# Add Files
set(LOCAL_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/is_bezier_spline.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/is_point.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/is_rational_bezier_spline.hpp
)

# Add All Subdirectories

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src/utils/type_traits)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)