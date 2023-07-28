set(SPLINEPY_SOURCES
    splinepy/bezier_explicit1.cpp
    splinepy/bezier_explicit2.cpp
    splinepy/bezier_explicit3.cpp
    splinepy/rational_bezier_explicit1.cpp
    splinepy/rational_bezier_explicit2.cpp
    splinepy/rational_bezier_explicit3.cpp
    splinepy/bspline_explicit1.cpp
    splinepy/bspline_explicit2.cpp
    splinepy/bspline_explicit3.cpp
    splinepy/nurbs_explicit1.cpp
    splinepy/nurbs_explicit2.cpp
    splinepy/nurbs_explicit3.cpp)
set(SPLINEPY_MORE_SOURCES
    splinepy/bezier_explicit4.cpp
    splinepy/bezier_explicit5.cpp
    splinepy/bezier_explicit6.cpp
    splinepy/bezier_explicit7.cpp
    splinepy/bezier_explicit8.cpp
    splinepy/bezier_explicit9.cpp
    splinepy/bezier_explicit10.cpp
    splinepy/rational_bezier_explicit4.cpp
    splinepy/rational_bezier_explicit5.cpp
    splinepy/rational_bezier_explicit6.cpp
    splinepy/rational_bezier_explicit7.cpp
    splinepy/rational_bezier_explicit8.cpp
    splinepy/rational_bezier_explicit9.cpp
    splinepy/rational_bezier_explicit10.cpp
    splinepy/bspline_explicit4.cpp
    splinepy/bspline_explicit5.cpp
    splinepy/bspline_explicit6.cpp
    splinepy/bspline_explicit7.cpp
    splinepy/bspline_explicit8.cpp
    splinepy/bspline_explicit9.cpp
    splinepy/bspline_explicit10.cpp
    splinepy/nurbs_explicit4.cpp
    splinepy/nurbs_explicit5.cpp
    splinepy/nurbs_explicit6.cpp
    splinepy/nurbs_explicit7.cpp
    splinepy/nurbs_explicit8.cpp
    splinepy/nurbs_explicit9.cpp
    splinepy/nurbs_explicit10.cpp)

set(EXPLICIT_SOURCES ${SPLINEPY_SOURCES})

if(SPLINEPY_MORE)
  set(EXPLICIT_SOURCES ${EXPLICIT_SOURCES} ${SPLINEPY_MORE_SOURCES})
endif(SPLINEPY_MORE)

if(SPLINEPY_BUILD_SHARED)
  add_library(explicit SHARED ${EXPLICIT_SOURCES})
else()
  add_library(explicit STATIC ${EXPLICIT_SOURCES})
endif(SPLINEPY_BUILD_SHARED)

# include
target_include_directories(
  explicit PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                   $<INSTALL_INTERFACE:${incl_dest}>)

# link
target_link_libraries(explicit PRIVATE bezman BSplineLib::splines napf)

# flags
target_compile_options(explicit PRIVATE ${SPLINEPY_FLAGS})

# defs
target_compile_definitions(explicit PRIVATE ${SPLINEPY_DEFS})

# features
target_compile_features(explicit PUBLIC cxx_std_17)

# install
install(
  TARGETS explicit
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION "${incl_dest}")
install(FILES "${project_config}" "${version_config}" DESTINATION "${cfg_dest}")
install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE "${namespace}"
  DESTINATION "${cfg_dest}")
