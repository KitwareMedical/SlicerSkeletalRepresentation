cmake_minimum_required(VERSION 3.13.4)

project(SkeletalRepresentation)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/KitwareMedical/SlicerSkeletalRepresentation/blob/master/README.md")
set(EXTENSION_CATEGORY "Shape Analysis")
set(EXTENSION_CONTRIBUTORS "Zhiyuan Liu (University of North Carolina at Chapel Hill), Junpyo Hong (University of North Carolina at Chapel Hill), Stephen M. Pizer (University of North Carolina at Chapel Hill), Jared Vicory (Kitware, Inc.), Pablo Hernandez-Cerdan (Kitware, Inc.), Beatriz Paniagua (Kitware, Inc.), Jean-Christophe Fillion-Robin (Kitware, Inc.), Connor Bowley (Kitware, Inc.)")
set(EXTENSION_DESCRIPTION "This extension handles skeletal representations of 3D data")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/KitwareMedical/SlicerSkeletalRepresentation/master/img/SRep.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/KitwareMedical/SlicerSkeletalRepresentation/master/img/VisualizationScreenshot.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

set(SUPERBUILD_TOPLEVEL_PROJECT inner)
# # Superbuild
# # Update CMake module path
set(EXTERNAL_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild)
set(EXTERNAL_PROJECT_FILE_PREFIX "External_")
set(LOCAL_PROJECT_NAME SkeletalRepresentation)
set(EXTENSION_NAME "${LOCAL_PROJECT_NAME}")


#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
mark_as_superbuild(Slicer_DIR)

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

#-----------------------------------------------------------------------------
# SuperBuild setup
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SRep)
add_subdirectory(SRepCreator)
add_subdirectory(SRepRefinement)
add_subdirectory(SRepWarper)
add_subdirectory(SRepPipelines)
add_subdirectory(SRepCreateAndRefine)
## NEXT_MODULE

#-----------------------------------------------------------------------------
set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
#list(APPEND EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS "${Foo_DIR};Foo;RuntimeLibraries;/")
set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

#-----------------------------------------------------------------------------
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
