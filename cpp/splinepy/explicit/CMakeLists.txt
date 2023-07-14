set(SPLINELIB_SOURCES
        splinelib/b_spline_explicit.cpp
        splinelib/nurbs_explicit.cpp
)

set(BEZMAN_SOURCES
        bezman/bezier_explicit.cpp
        bezman/rational_bezier_explicit.cpp
)

set(SPLINEPY_SOURCES
        splinepy/bezier_explicit.cpp
        splinepy/rational_bezier_explicit.cpp
        splinepy/bspline_explicit.cpp
        splinepy/nurbs_explicit.cpp
)

set(EXPLICIT_SOURCES
        ${SPLINELIB_SOURCES}
        ${BEZMAN_SOURCES}
        ${SPLINEPY_SOURCES}
)

if(SPLINEPY_MORE)
    set(EXPLICIT_SOURCES
            ${EXPLICIT_SOURCES}
            ${SPLINELIB_MORE_SOURCES}
            ${BEZMAN_MORE_SOURCES}
            ${SPLINEPY_MORE_SOURCES}
    )
endif(SPLINEPY_MORE)

if(SPLINEPY_BUILD_SHARED)
    add_library(explicit SHARED ${EXPLICIT_SOURCES})
else()
    add_library(explicit STATIC ${EXPLICIT_SOURCES})
endif(SPLINEPY_BUILD_SHARED)

# include
target_include_directories(explicit
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>
        $<INSTALL_INTERFACE:${incl_dest}>)

# link
target_link_libraries(explicit
    PRIVATE
        bezman
        BSplineLib::splines
        napf)

# flags
target_compile_options(explicit
    PRIVATE
        ${SPLINEPY_FLAGS})

# defs
target_compile_definitions(explicit
    PRIVATE
        ${SPLINEPY_DEFS})

# features
target_compile_features(explicit
    PUBLIC
        cxx_std_17)

# install
install(
    TARGETS explicit
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
