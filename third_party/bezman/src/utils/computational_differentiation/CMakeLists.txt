# Add Files
set(LOCAL_HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/algo_diff_type.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/algo_diff_type.inc
)

# Add All Subdirectories

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src/utils/computational_differentiation)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)