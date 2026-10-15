# Paths to the real DTDs (choose whatever is correct for your tree)
set(CONVEX_DTD "Convex.dtd")
set(CONTACT_DTD "Contact.dtd")


# Generate the configured DTD into the build dir (or wherever you want)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Contact.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Contact.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Convex.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Convex.dtd"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains3D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains3D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains3D.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains3D.dtd"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Generic_Grains3D.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Generic_Grains3D.dtd"
  @ONLY
)
