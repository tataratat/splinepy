add_subdirectory(fitting)

if(SPLINEPY_COMPILE_PYTHON)
    add_subdirectory(py)
else()
    # make Py<spline-class> available
    add_library(splinepy INTERFACE)
    add_library(splinepy::splinepy ALIAS splinepy)

    # include
    target_include_directories(splinepy
        INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>)
            #$<INSTALL_INTERFACE:${incl_dest}>
 
    # link
    target_link_libraries(splinepy
        INTERFACE
            pybind11::headers
            fitting
            bezman
            napf
            SplineLib::splines
            SplineLib::input_output)
endif()
