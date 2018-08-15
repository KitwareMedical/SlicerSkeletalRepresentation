cmake_minimum_required(VERSION 3.5)

project(SkeletalRepresentation)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SkeletalRepresentation")
set(EXTENSION_CATEGORY "Topology, Skeletonization")
set(EXTENSION_CONTRIBUTORS "Jared Vicory (Kitware, Inc.), Pablo Hernandez-Cerdan (Kitware, Inc.), Beatriz Paniagua (Kitware, Inc.)")
set(EXTENSION_DESCRIPTION "Collection of Skeletal Model Extensions")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SkeletalRepresentation.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SkeletalRepresentation/Screenshots/1.png")
set(EXTENSION_DEPENDS Eigen3) # Specified as a space separated string, a list or 'NA' if any

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
