set(SRCS_FITTING fitting.cpp)

add_library(fitting OBJECT ${SRCS_FITTING})
target_include_directories(fitting PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(fitting PRIVATE ${SPLINELIBPY_FLAGS})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(fitting PRIVATE ${GNU_FLAGS_R})
endif()
