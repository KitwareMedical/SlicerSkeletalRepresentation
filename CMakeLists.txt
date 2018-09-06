cmake_minimum_required(VERSION 3.5)

project(SkeletalRepresentation)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SkeletalRepresentation")
set(EXTENSION_CATEGORY "Topology, Skeletonization")
set(EXTENSION_CONTRIBUTORS "Zhiyuan Liu (University of North Carolina at Chapel Hill), Junpyo Hong (University of North Carolina at Chapel Hill), Stephen M. Pizer (University of North Carolina at Chapel Hill), Jared Vicory (Kitware, Inc.), Pablo Hernandez-Cerdan (Kitware, Inc.), Beatriz Paniagua (Kitware, Inc.), Jean-Christophe Fillion-Robin (Kitware, Inc.)")
set(EXTENSION_DESCRIPTION "Toolkit that provides users to establish, refine, and visualize skeletal representation of an object.")
set(EXTENSION_ICONURL "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SkeletalRepresentation.png")
set(EXTENSION_SCREENSHOTURLS "https://github.com/ZhiyLiu/SRepExtension/blob/master/visualize_s_rep.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SkeletalRepresentationVisualizer)
add_subdirectory(SkeletalRepresentationInitializer)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
