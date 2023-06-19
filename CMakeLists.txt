cmake_minimum_required(VERSION 3.12.0)

project(splinepy VERSION 0.0.2 LANGUAGES CXX)

# build options
option(SPLINEPY_COMPILE_SPLINELIB "Compile splinelib together." ON)
option(SPLINEPY_VERBOSE_MAKE
        "Verbose `make` output. Alias to CMAKE_VERBOSE_MAKEFILE" OFF)
option(SPLINEPY_MORE "Compile a full set of splines." OFF)
option(SPLINEPY_BUILD_SHARED "build shared library for splinepy" OFF)
option(SPLINEPY_COMPILE_PYTHON "Compile python module." ON)
option(SPLINEPY_ENABLE_WARNINGS "Add warning flags" OFF)
option(SPLINEPY_BUILD_EXPLICIT "Explicit instantiation of (mainly third party) template classes" ON)


# config
set(exe_dest "bin")
set(incl_dest "include")
set(lib_dest "lib")
set(cfg_dest "${lib_dest}/cmake/${PROJECT_NAME}")
set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(version_config "${gen_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${gen_dir}/${PROJECT_NAME}Config.cmake")
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(namespace "${PROJECT_NAME}::")


# global flags
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(SPLINEPY_DEFS $<$<NOT:$<CONFIG:Debug>>:NDEBUG>)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(SPLINEPY_FLAGS -fPIC)
  set(SPLINEPY_OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(SPLINEPY_OPTIMIZATION_FLAGS
        $<$<NOT:$<CONFIG:Debug>>:/O2>
        $<$<NOT:$<CONFIG:Debug>>:/w>)
  set(SPLINEPY_FLAGS /bigobj)
endif()

# compiler specific flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(SPLINEPY_WARNING_FLAGS
        -Wall -Wextra -Wpedantic -Wzero-as-null-pointer-constant -Wno-unused)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(SPLINEPY_WARNING_FLAGS
        -Wall -Wextra -Wmost -Wextra -Wpedantic -Wunreachable-code
        -Wshadow -Wfloat-equal -Weffc++ -Wno-unused-parameter
        -Wno-unused-variable -Wzero-as-null-pointer-constant)
else()
  message(WARNING
            "Unsupported compiler."
            "splinepy is tested with GNU and Clang.")
endif()

if(SPLINEPY_ENABLE_WARNINGS)
  set(SPLINEPY_FLAGS ${SPLINEPY_FLAGS} ${SPLINEPY_WARNING_FLAGS})
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(SPLINEPY_FLAGS ${SPLINEPY_FLAGS} ${SPLINEPY_OPTIMIZATION_FLAGS})
endif()

if(SPLINEPY_VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if(SPLINEPY_MORE)
  set(SPLINEPY_DEFS ${SPLINEPY_DEFS} SPLINEPY_MORE)
endif(SPLINEPY_MORE)

if(SPLINEPY_BUILD_EXPLICIT)
    set(SPLINEPY_DEFS ${SPLINEPY_DEFS}
            SPLINEPY_BUILD_EXPLICIT)
endif(SPLINEPY_BUILD_EXPLICIT)

add_subdirectory(third_party)
add_subdirectory(cpp/splinepy)
