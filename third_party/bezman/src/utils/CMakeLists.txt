# Add Files
set(LOCAL_HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/base64.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/binomialcoefficientlookuptablecreator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/export.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/export.inc
  ${CMAKE_CURRENT_SOURCE_DIR}/fastbinomialcoefficient.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/logger.hpp
)

# Add All Subdirectories
add_subdirectory(computational_differentiation)
add_subdirectory(type_traits)
add_subdirectory(algorithms)

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src/utils)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)