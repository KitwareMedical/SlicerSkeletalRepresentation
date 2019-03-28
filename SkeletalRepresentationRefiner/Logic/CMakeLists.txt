project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  vtkSlicerSkeletalRepresentationInterpolater.h
  vtkSlicerSkeletalRepresentationInterpolater.cpp
  vtkSpoke.h
  vtkSpoke.cpp
  vtkSrep.h
  vtkSrep.cpp
  newuoa.h
  vtkPolyData2ImageData.cpp
  vtkPolyData2ImageData.h
  vtkImage2SignedDistanceMap.h
  vtkImage2SignedDistanceMap.cpp
  vtkAntiAlias.h
  vtkAntiAlias.cpp
  vtkApproximateSignedDistanceMap.h
  vtkApproximateSignedDistanceMap.cpp
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
