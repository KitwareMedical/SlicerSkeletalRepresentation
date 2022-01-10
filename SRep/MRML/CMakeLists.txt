project(vtkSlicer${MODULE_NAME}ModuleMRML)

find_package(RapidJSON REQUIRED)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_MRML_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  ${RapidJSON_INCLUDE_DIR}
  ${vtkSlicerSRepModuleMRML_INCLUDE_DIRS}
  )

set(${KIT}_SRCS
  vtkMRMLEllipticalSRepNode.cxx
  vtkMRMLRectangularGridSRepNode.cxx
  vtkMRMLSRepNode.cxx
  vtkMRMLSRepDisplayNode.cxx
  vtkMRMLSRepStorageNode.cxx

  # non MRML nodes
  vtkEllipticalSRep.cxx
  vtkEllipticalSRep.h
  vtkMeshSRepInterface.h
  vtkSRepSkeletalPoint.cxx
  vtkSRepSkeletalPoint.h
  vtkSRepSpoke.cxx
  vtkSRepSpoke.h
  vtkSRepSpokeMesh.cxx
  vtkSRepSpokeMesh.h
  )

set(${KIT}_TARGET_LIBRARIES
  ${MRML_LIBRARIES}
  srep::srep
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleMRML(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
