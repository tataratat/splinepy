# build fitting first - currently very independent
add_subdirectory(fitting)


## create splinepy target - enables cpp standalone use
# define srcs
set(SPLINEPY_SRCS
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier1.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier2.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier3.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier1.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier2.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier3.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline1.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline2.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline3.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs1.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs2.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs3.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/splinepy_base.cpp
)

set(SPLINEPY_MORE_SRCS
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier4.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier5.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier6.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier7.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier8.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier9.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bezier10.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier4.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier5.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier6.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier7.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier8.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier9.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/rational_bezier10.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline4.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline5.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline6.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline7.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline8.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline9.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/bspline10.cpp

        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs4.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs5.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs6.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs7.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs8.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs9.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create/nurbs10.cpp
)

if(SPLINEPY_MORE)
  set(SPLINEPY_SRCS ${SPLINEPY_SRCS} ${SPLINEPY_MORE_SRCS})
endif(SPLINEPY_MORE)


# target
if(SPLINEPY_BUILD_SHARED)
    add_library(splinepy SHARED ${SPLINEPY_SRCS})
else()
    add_library(splinepy STATIC ${SPLINEPY_SRCS})
endif(SPLINEPY_BUILD_SHARED)
# alias for convenience
add_library(splinepy::splinepy ALIAS splinepy)


# include
target_include_directories(splinepy
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>
        $<INSTALL_INTERFACE:${incl_dest}>)


# link
target_link_libraries(splinepy
    PUBLIC
        fitting
        bezman
        napf
        SplineLib::splines
        SplineLib::input_output)

# explicit?
if (SPLINEPY_BUILD_EXPLICIT)
add_subdirectory(explicit)
target_link_libraries(splinepy
    PUBLIC
        explicit)
endif()


# flags
target_compile_options(splinepy
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

# python?
if(SPLINEPY_COMPILE_PYTHON)
    add_subdirectory(py)
endif()


# cmake config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${cfg_dest}"
)


# install
install(
    TARGETS splinepy
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
# .hpp files - keep tree
install(DIRECTORY "${PROJECT_SOURCE_DIR}/cpp/splinepy"
        DESTINATION ${incl_dest} # target directory
        FILES_MATCHING PATTERN "*.hpp"
)
