set(SRCS_FITTING fitting.cpp)

if(SPLINEPY_BUILD_SHARED)
    add_library(fitting SHARED ${SRCS_FITTING})
else()
    add_library(fitting STATIC ${SRCS_FITTING})
endif()
target_include_directories(fitting PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# 11 is enough, but 17 allows nested namespace..
target_compile_features(fitting PUBLIC cxx_std_17)

target_compile_options(fitting PRIVATE ${SPLINEPY_FLAGS})
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(fitting PRIVATE ${SPLINEPY_OPTIMIZAION_FLAGS})
endif()

install(
    TARGETS fitting
    EXPORT ${TARGETS_EXPORT_NAME}
    LIBRARY DESTINATION ${lib_dest}
    LIBRARY DESTINATION ${lib_dest}
    INCLUDES DESTINATION ${incl_dest}
)
