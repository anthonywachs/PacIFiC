cmake_minimum_required(VERSION 3.12)

find_program(PACIFIC_DOXYGEN_EXECUTABLE NAMES doxygen)
find_program(PACIFIC_MKDOCS_EXECUTABLE NAMES mkdocs)

if(NOT PACIFIC_DOXYGEN_EXECUTABLE)
  message(WARNING "doxygen not found: docs_doxygen/docs_build/docs_develop targets will fail.")
endif()

if(NOT PACIFIC_MKDOCS_EXECUTABLE)
  message(WARNING "mkdocs not found: docs_build/docs_develop targets will fail.")
endif()

set(PACIFIC_DOCS_WORKDIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(PACIFIC_DOCS_DOXYFILES
  DoxyFLUID
  DoxyMAC
  DoxyGrains3D
)

add_custom_target(docs_doxygen
  COMMAND "${PACIFIC_DOXYGEN_EXECUTABLE}" "${PACIFIC_DOCS_WORKDIR}/DoxyFLUID"
  COMMAND "${PACIFIC_DOXYGEN_EXECUTABLE}" "${PACIFIC_DOCS_WORKDIR}/DoxyMAC"
  COMMAND "${PACIFIC_DOXYGEN_EXECUTABLE}" "${PACIFIC_DOCS_WORKDIR}/DoxyGrains3D"
  WORKING_DIRECTORY "${PACIFIC_DOCS_WORKDIR}"
  COMMENT "Generating Doxygen documentation"
  VERBATIM
  USES_TERMINAL
)

add_custom_target(docs_clean
  COMMAND "${CMAKE_COMMAND}" -E rm -rf
    "${PACIFIC_DOCS_WORKDIR}/site"
    "${PACIFIC_DOCS_WORKDIR}/pages/doxygen/grains3d_html"
    "${PACIFIC_DOCS_WORKDIR}/pages/doxygen/fluid_html"
    "${PACIFIC_DOCS_WORKDIR}/pages/doxygen/mac_html"
  WORKING_DIRECTORY "${PACIFIC_DOCS_WORKDIR}"
  COMMENT "Cleaning generated docs outputs"
  VERBATIM
)

add_custom_target(docs_build
  COMMAND "${PACIFIC_MKDOCS_EXECUTABLE}" build
  WORKING_DIRECTORY "${PACIFIC_DOCS_WORKDIR}"
  COMMENT "Building MkDocs site"
  VERBATIM
  USES_TERMINAL
)

add_dependencies(docs_build docs_clean docs_doxygen)

add_custom_target(docs_develop
  COMMAND "${PACIFIC_MKDOCS_EXECUTABLE}" serve
  WORKING_DIRECTORY "${PACIFIC_DOCS_WORKDIR}"
  COMMENT "Serving MkDocs site locally"
  VERBATIM
  USES_TERMINAL
)

add_dependencies(docs_develop docs_doxygen)
