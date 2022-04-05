project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  SRepInterpolater.h
  SRepInterpolater.cpp
  Spoke.h
  Spoke.cpp
  SRep.h
  SRep.cpp
  newuoa.h
  PolyData2ImageData.cpp
  PolyData2ImageData.h
  ApproximateSignedDistanceMap.h
  ApproximateSignedDistanceMap.cpp
  GradientDistanceFilter.cpp
  GradientDistanceFilter.h
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  VTK::eigen
  vtkSlicerMarkupsModuleMRML
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
