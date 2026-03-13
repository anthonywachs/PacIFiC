function(octree_add_executable SOURCE_FILE)
  get_filename_component(source_name ${SOURCE_FILE} NAME_WE)
  set(output_c "${CMAKE_BINARY_DIR}/_${source_name}.c")

  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/_${source_name}.c"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SOURCE_FILE}" "${CMAKE_BINARY_DIR}/${source_name}.c"
    COMMAND $<TARGET_FILE:basilisk::qcc>
      "${source_name}.c"
      -I"${CMAKE_SOURCE_DIR}"
      -I"${Octree_DLMFD_INCLUDE_DIR}"
      -I"${Octree_VTKHYPERTREE_INCLUDE_DIR}"
      -DTRACE=2
      -source
    DEPENDS ${SOURCE_FILE}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" 
  )

  add_executable(${source_name} "_${source_name}.c")

  target_include_directories(${source_name} PRIVATE "${CMAKE_SOURCE_DIR}")

  target_link_libraries(${source_name}
    PRIVATE
      PacIFiC::Grains3D
      MPI::MPI_C
      HDF5::HDF5
      m
  )

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )

endfunction()

