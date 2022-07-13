set(SRCS_FITTING fitting.cpp)

add_library(fitting OBJECT ${SRCS_FITTING})
target_include_directories(fitting PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(fitting PRIVATE ${SPLINEPY_FLAGS})
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(fitting PRIVATE ${SPLINEPY_OPTIMIZAION_FLAGS})
endif()
