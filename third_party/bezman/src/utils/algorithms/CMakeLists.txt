# Add Files
set(LOCAL_HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/hypercube.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/int_power.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sort.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/point_uniquifier.hpp
)

# Add All Subdirectories

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src/utils/algorithms)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)