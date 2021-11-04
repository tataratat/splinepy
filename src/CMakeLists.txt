set(SRCS splinelib.cpp)

pybind11_add_module(splinelibpy ${SRCS})
target_include_directories(splinelibpy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(splinelibpy PRIVATE ${SPLINELIBPY_FLAGS})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(splinelibpy PRIVATE ${GNU_FLAGS_R})
endif()

# Link SplineLib
if(COMPILE_SPLINELIB)
    target_link_libraries(splinelibpy PRIVATE splines input_output)
else()
    find_package(SplineLib REQUIRED)
    target_link_libraries(splinelibpy 
        PRIVATE
            SplineLib::splines
            SplineLib::input_output
    )
endif()

# Link Fitting
add_subdirectory(fitting)
target_link_libraries(splinelibpy PRIVATE fitting)
