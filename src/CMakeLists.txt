# create splinepy target - enables cpp standalone use define srcs
set(SPLINEPY_SRCS
    ${PROJECT_SOURCE_DIR}/src/proximity/proximity.cpp
    ${PROJECT_SOURCE_DIR}/src/utils/coordinate_pointers.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/helpers/extract.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier1.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier2.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier3.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier1.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier2.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier3.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/splinepy_base.cpp)

set(SPLINEPY_MORE_SRCS
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier4.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier5.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier6.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier7.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier8.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier9.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/bezier10.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier4.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier5.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier6.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier7.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier8.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier9.cpp
    ${PROJECT_SOURCE_DIR}/src/splines/create/rational_bezier10.cpp)

# extend src for high dim splines
if(SPLINEPY_MORE)
  set(SPLINEPY_SRCS ${SPLINEPY_SRCS} ${SPLINEPY_MORE_SRCS})
endif(SPLINEPY_MORE)

# target
add_library(splinepy ${SPLINEPY_LIB_TYPE} ${SPLINEPY_SRCS})
# alias for convenience
add_library(splinepy::splinepy ALIAS splinepy)

# include
target_include_directories(
  splinepy PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                  $<INSTALL_INTERFACE:${incl_dest}>)

# make sure BSplineLib is available
if(NOT TARGET BSplineLib::splines OR NOT TARGET BSplineLib::input_output)
  find_package(BSplineLib REQUIRED)
endif()

# link
target_link_libraries(splinepy PUBLIC bezman napf BSplineLib::splines
                                      BSplineLib::input_output)

# explicit?
if(SPLINEPY_BUILD_EXPLICIT)
  add_subdirectory(explicit)
  target_link_libraries(splinepy PUBLIC explicit)
endif()

# flags
target_compile_options(splinepy PRIVATE ${SPLINEPY_FLAGS})

# defs
target_compile_definitions(splinepy PRIVATE ${SPLINEPY_DEFS})

# features
target_compile_features(splinepy PUBLIC cxx_std_17)

# python?
if(SPLINEPY_COMPILE_PYTHON)
  add_subdirectory(py)
endif()

# cmake config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${version_config}"
                                 COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in" "${project_config}"
  INSTALL_DESTINATION "${cfg_dest}")

# install main target
install(
  TARGETS splinepy
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION "${incl_dest}")

# config files
install(FILES "${project_config}" "${version_config}" DESTINATION "${cfg_dest}")
install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE "${namespace}"
  DESTINATION "${cfg_dest}")

# header files - keep tree
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/include/splinepy"
  DESTINATION ${incl_dest} # target directory
  FILES_MATCHING
  PATTERN "*.hpp")

# inline files - keep tree
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/include/splinepy"
  DESTINATION ${incl_dest} # target directory
  FILES_MATCHING
  PATTERN "*.inl")
