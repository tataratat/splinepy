set(XNS_HELPER_C_SRCS ${PROJECT_SOURCE_DIR}/src/xns/splinepy_scene.cpp)

set(XNS_HELPER_F_SRCS ${PROJECT_SOURCE_DIR}/src/xns/splinepy_scene.f90)

# C/CPP
add_library(xns_helper_c ${SPLINEPY_LIB_TYPE} ${XNS_HELPER_C_SRCS})
add_library(splinepy::xns_helper_c ALIAS xns_helper_c)
target_link_libraries(xns_helper_c PUBLIC splinepy::splinepy_python)

# F
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/splinepy_f_modules)
add_library(xns_helper_f ${SPLINEPY_LIB_TYPE} ${XNS_HELPER_F_SRCS})
target_link_libraries(xns_helper_f PUBLIC xns_helper_c)

install(
  TARGETS xns_helper_c
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION ${include_dest})
install(
  TARGETS xns_helper_f
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION ${include_dest})
file(GLOB_RECURSE SPLINEPY_F_MOD "${CMAKE_Fortran_MODULE_DIRECTORY}/*.mod")
install(FILES ${SPLINEPY_F_MOD} DESTINATION ${include_dest})
install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE ${namespace}
  DESTINATION ${cfg_dest})
