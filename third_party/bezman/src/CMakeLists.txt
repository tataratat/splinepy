# Add files
set(LOCAL_HEADERS 
      ${CMAKE_CURRENT_SOURCE_DIR}/bezier_spline.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/bezier_spline.inc
      ${CMAKE_CURRENT_SOURCE_DIR}/bezier_spline_group.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/bezier_spline_group.inc
      ${CMAKE_CURRENT_SOURCE_DIR}/point.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/rational_bezier_spline.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/rational_bezier_spline.inc
)

# Add All Subdirectories
add_subdirectory(utils)

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)