
# Find all *.exec in this source directory
file(GLOB EXEC_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.exec")

# Copy them into the matching build directory
foreach(src ${EXEC_SCRIPTS})
  get_filename_component(name "${src}" NAME)
  configure_file("${src}" "${CMAKE_CURRENT_BINARY_DIR}/${name}" COPYONLY)
endforeach()