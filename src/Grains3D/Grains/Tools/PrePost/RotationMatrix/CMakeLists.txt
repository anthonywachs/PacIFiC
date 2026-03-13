add_executable(Grains3D_rotationmatrix
  src/PP_rotationmatrix.cpp
)

set_target_properties(Grains3D_rotationmatrix PROPERTIES OUTPUT_NAME rotationmatrix)

target_link_libraries(Grains3D_rotationmatrix
  PRIVATE
    Grains3D
)

