cmake_minimum_required(VERSION 3.13.4)

project(SRep2)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SRep2")
set(EXTENSION_CATEGORY "SRep")
set(EXTENSION_CONTRIBUTORS "Connor Bowley (Kitware, Inc.)")
set(EXTENSION_DESCRIPTION "This extension handles SRep representations of data")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SRep2.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SRep2/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

add_subdirectory(CommonLibrary)

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SRep)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
