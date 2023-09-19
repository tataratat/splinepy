set(XNS_HELPER_F_SRCS ${PROJECT_SOURCE_DIR}/src/xns/splinepy_scene.f90)

# F - it's just an interface, so .so itself will be an empty one
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/splinepy_f_modules)
add_library(xns_helper_f ${SPLINEPY_LIB_TYPE} ${XNS_HELPER_F_SRCS})
target_link_libraries(xns_helper_f PUBLIC splinepy_python)

install(
  TARGETS xns_helper_f
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${lib_dest}
  ARCHIVE DESTINATION ${lib_dest}
  INCLUDES
  DESTINATION ${incl_dest})
file(GLOB_RECURSE SPLINEPY_F_MOD "${CMAKE_Fortran_MODULE_DIRECTORY}/*.mod")
install(FILES ${SPLINEPY_F_MOD} DESTINATION ${incl_dest})
install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE ${namespace}
  DESTINATION ${cfg_dest})
