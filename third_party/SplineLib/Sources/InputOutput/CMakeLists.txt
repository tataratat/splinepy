# Copyright (c) 2018–2021 SplineLib
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS
# OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pkg_search_module(PUGIXML REQUIRED QUIET pugixml)

set(
    HEADERS
    iges.hpp
    irit.hpp
    operations.hpp
    operations.inc
    vtk.hpp
    xml.hpp
)

set(
    SOURCES
    iges.cpp
    operations.cpp
    irit.cpp
    vtk.cpp
    xml.cpp
#
    ${HEADERS}
)
set_source_files_properties(${HEADERS} PROPERTIES LANGUAGE CXX HEADER_FILE_ONLY TRUE)

add_library(input_output ${SOURCES})
target_include_directories(input_output PUBLIC ${INCLUDE_DIRECTORIES})
target_include_directories(input_output SYSTEM PRIVATE ${PUGIXML_INCLUDEDIR})
target_link_libraries(input_output PUBLIC utilities splines 
                                   PRIVATE ${PUGIXML_LDFLAGS} "-Wl,-rpath,${PUGIXML_LIBRARY_DIRS}")
target_compile_definitions(input_output INTERFACE ${COMPILE_DEFINITIONS} PRIVATE ${PUGIXML_CFLAGS_OTHER})
target_compile_options(input_output PRIVATE ${COMPILE_OPTIONS})

install(
        TARGETS input_output
        EXPORT ${targets_export_name}
        LIBRARY DESTINATION ${library_install_directory}
        ARCHIVE DESTINATION ${library_install_directory}/static
        RUNTIME DESTINATION ${executable_install_directory}
        INCLUDES DESTINATION ${include_install_directory}
)

install(FILES ${HEADERS} DESTINATION ${include_install_directory_headers}/Sources/InputOutput)
