add_library(srep
  EllipticalSRep.cpp
  Point3d.cpp
  ReadHeaderFile.cpp
  RectangularGridSRep.cpp
  SkeletalPoint.cpp
  Spoke.cpp
  SpokeMesh.cpp
  SRepIO.cpp
  Vector3d.cpp
)
target_include_directories(srep PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(srep
PUBLIC
  -lexpat
  ${VTK_LIBRARIES}
)
add_library(srep::srep ALIAS srep)

#TODO: Add testing. Need to put GTest in superbuild for this
# if (${BUILD_TESTING})
#   add_subdirectory(testing)
# endif()
