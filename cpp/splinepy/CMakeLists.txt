# build fitting first - currently very independent
# also only used py python side
add_subdirectory(fitting)

## create splinepy target - enables cpp standalone use
# define srcs
set(SPLINEPY_SRCS
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create_bezier.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create_rational_bezier.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create_bspline.cpp
        ${PROJECT_SOURCE_DIR}/cpp/splinepy/splines/create_nurbs.cpp
)

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
        #pybind11::headers  # TODO is this necessary?
        fitting
        bezman
        napf
        SplineLib::splines
        SplineLib::input_output)

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
