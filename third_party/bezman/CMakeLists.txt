####
# General information
####
# Required CMAKE 
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Define Project
project("bezman"
        VERSION 0.0.1
        DESCRIPTION "A header library to facilitate the manipulation of Bezier-splines"
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(namespace "${PROJECT_NAME}::")
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(MY_FANCY_TARGET_NAME "${PROJECT_NAME}")

####
# Global Options
####
# Testing
option(GOOGLETEST "Build with GoogleTest" OFF)

# Allow Logging
option(LOGGING "Enable Logging" OFF)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Print out information about compiler
message(STATUS "Compile using ${CMAKE_CXX_COMPILER_ID} compiler")

# Set compile information
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3 -march=native>)
    set(RUNTIME_CHECKS_DEBUG -D_GLIBCXX_DEBUG)
    set(WARNING_FLAGS -Wall -Wextra -Wpedantic -Wzero-as-null-pointer-constant $<$<NOT:$<CONFIG:Debug>>:-Wno-unused>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("CLANG generally not supported, and not sufficiently tested")
    set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3 -march=native>)
    set(RUNTIME_CHECKS_DEBUG)
    set(WARNING_FLAGS -Wall -Wextra -Wmost -Wextra -Wpedantic -Wunreachable-code -Wshadow -Wfloat-equal -Weffc++
                       $<$<NOT:$<CONFIG:Debug>>:-Wno-unused-parameter -Wno-unused-variable>
                       $<$<BOOL:${GOOGLETEST}>:-Wno-gnu-zero-variadic-macro-arguments>)
else()
    message(ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

####
# Add contents
####

# Aux Variables (files will be directly set to this location)
set(include_destination "include")
set(include_destination_headers ${include_destination}/${PROJECT_NAME})

# Add Files by going through directory structure and add support for IDEs
add_subdirectory(src)
set_source_files_properties(${HEADERS} PROPERTIES LANGUAGE CXX HEADER_FILE_ONLY TRUE)
source_group(TREE ${PROJECT_SOURCE_DIR} FILES ${HEADERS})

####
# Packaging
####
set(lib_destination "lib")
set(config_destination "${lib_destination}/cmake/${PROJECT_NAME}")
set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(version_config "${gen_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${gen_dir}/${PROJECT_NAME}Config.cmake")

# Load cmake Macros
include(CMakePackageConfigHelpers)

# Set Up Version
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion
    )

# Automatic configuration file
configure_package_config_file(
    "cmake/config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_destination}"
    )

# Define Target
add_library(${MY_FANCY_TARGET_NAME} INTERFACE)
target_include_directories(${MY_FANCY_TARGET_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/..>)

## INSTALLATION PROCESS
install(TARGETS ${MY_FANCY_TARGET_NAME} EXPORT ${TARGETS_EXPORT_NAME} INCLUDES DESTINATION ${include_destination})
install(FILES ${project_config} ${version_config} DESTINATION ${config_destination})
install(EXPORT ${TARGETS_EXPORT_NAME} NAMESPACE ${namespace} DESTINATION ${config_destination})

### User Output
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Build debug!")
    target_compile_definitions(${MY_FANCY_TARGET_NAME} INTERFACE NDEBUG)
else()
    message(STATUS "Build release!")
endif()

# Install and provide testing 
if(GOOGLETEST)
    message(STATUS "Build tests!")
    
    # Allow Logging
    set(LOGGING OFF)
    

    # The following 10 lines are from the google docs documentation to ensure version compatibility
    # Fetch the newest gtest version and make it locally available
    include(FetchContent)
    FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
    )
    
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(tests)
else()
    message(STATUS "No tests!")
endif()

# Enable Logging by setting a preprocessor flag
if (LOGGING)
    target_compile_definitions(${MY_FANCY_TARGET_NAME} INTERFACE ENABLE_LOGGING)
endif()
