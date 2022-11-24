set(SPLINELIB_SOURCES
        splinelib/b_spline_explicit.cpp
        splinelib/nurbs_explicit.cpp
)

set(BEZMAN_SOURCES
        bezman/bezier_explicit.cpp
        bezman/rational_bezier_explicit.cpp
)

set(EXPLICIT_SOURCES
        ${SPLINELIB_SOURCES}
#        ${BEZMAN_SOURCES}
)

if(SPLINEPY_BUILD_SHARED)
    add_library(third_party_explicit SHARED ${EXPLICIT_SOURCES})
else()
    add_library(third_party_explicit STATIC ${EXPLICIT_SOURCES})
endif(SPLINEPY_BUILD_SHARED)

# include
target_include_directories(third_party_explicit
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>
        $<INSTALL_INTERFACE:${incl_dest}>)

# link
target_link_libraries(third_party_explicit
    PRIVATE
        bezman
        SplineLib::splines)

# flags
target_compile_options(third_party_explicit
    PRIVATE
        ${SPLINEPY_FLAGS})

# defs
target_compile_definitions(splinepy
    PRIVATE
        ${SPLINEPY_DEFS})

# features
target_compile_features(splinepy
    PUBLIC
        cxx_std_17)

# install
install(
    TARGETS third_party_explicit
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION ${lib_dest}
    ARCHIVE DESTINATION ${lib_dest}
    INCLUDES DESTINATION "${incl_dest}"
)
install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${cfg_dest}"
)
install(
    EXPORT "${TARGETS_EXPORT_NAME}"
    NAMESPACE "${namespace}"
    DESTINATION "${cfg_dest}"
)
