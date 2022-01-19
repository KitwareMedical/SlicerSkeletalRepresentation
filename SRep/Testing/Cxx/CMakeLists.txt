include(GoogleTest)
message("CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")
message("GTest_DIR ${GTest_DIR}")

find_package(GTest REQUIRED CONFIG)

add_executable(qSlicerSRepModuleUnitTests
  Point3dTest.cxx
  SkeletalPointTest.cxx
  SpokeTest.cxx
  Vector3dTest.cxx
)

target_link_libraries(qSlicerSRepModuleUnitTests
  vtkSlicerSRepModuleMRML
  GTest::gtest_main
)

add_test(NAME qSlicerSRepModuleUnitTests COMMAND qSlicerSRepModuleUnitTests)
set_property(TEST qSlicerSRepModuleUnitTests PROPERTY LABELS qSlicerSRepModule)
