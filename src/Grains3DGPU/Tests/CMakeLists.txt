set(GRAINS3D_GPU_TEST_CANDIDATES
  geometry/test_vector3.cpp
  geometry/test_matrix3.cpp
  geometry/test_quaternion.cpp
  geometry/RotationMathTest.cpp
  geometry/test_transform3.cpp
  geometry/test_transform.cpp
  collision/test_gjk.cpp
  collision/test_obb.cpp
  collision/test_obc.cpp
  collision/test_collision_detection.cpp
  collision/test_neighborlist.cpp
  component/test_rigid_body.cpp
  component/test_convex.cpp
  contactForce/test_contact_table_host.cpp
  contactForce/test_contact_table_device.cpp
)

set(GRAINS3D_GPU_TEST_SRCS main.cpp)
foreach(test_src IN LISTS GRAINS3D_GPU_TEST_CANDIDATES)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_src}")
    list(APPEND GRAINS3D_GPU_TEST_SRCS "${test_src}")
  endif()
endforeach()

add_executable(Grains3DGPU_tests ${GRAINS3D_GPU_TEST_SRCS})

set_source_files_properties(${GRAINS3D_GPU_TEST_SRCS} PROPERTIES LANGUAGE CUDA)

set_target_properties(Grains3DGPU_tests PROPERTIES
  CUDA_STANDARD ${GRAINS3D_GPU_CUDA_STANDARD}
  CUDA_STANDARD_REQUIRED ${GRAINS3D_GPU_CUDA_STANDARD_REQUIRED}
  CUDA_ARCHITECTURES ${GRAINS3D_GPU_CUDA_ARCHITECTURES}
  CUDA_SEPARABLE_COMPILATION ${GRAINS3D_GPU_CUDA_SEPARABLE_COMPILATION}
  CUDA_RESOLVE_DEVICE_SYMBOLS ${GRAINS3D_GPU_CUDA_RESOLVE_DEVICE_SYMBOLS}
)

foreach(cuda_opt IN LISTS GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST)
  target_compile_options(Grains3DGPU_tests PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:${cuda_opt}>
  )
endforeach()

target_compile_definitions(Grains3DGPU_tests PRIVATE TEST_BUILD=1)

find_package(GTest REQUIRED)

set(_GTEST_LIBS "")
if(TARGET GTest::gtest AND TARGET GTest::gtest_main)
  list(APPEND _GTEST_LIBS GTest::gtest GTest::gtest_main)
elseif(TARGET GTest::GTest AND TARGET GTest::Main)
  list(APPEND _GTEST_LIBS GTest::GTest GTest::Main)
else()
  message(FATAL_ERROR "Unsupported GTest target names from find_package(GTest)")
endif()

target_link_libraries(Grains3DGPU_tests PRIVATE
  Grains3DGPU_obj
  ${_GTEST_LIBS}
)

include(GoogleTest)
gtest_discover_tests(Grains3DGPU_tests
  TEST_PREFIX "Grains3DGPU."
  DISCOVERY_TIMEOUT 60
)
