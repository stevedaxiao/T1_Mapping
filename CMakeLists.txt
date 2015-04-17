cmake_minimum_required(VERSION 2.8.9)

project(T1_Mapping)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://wiki.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/T1_Mapping")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Xiao Da (MGH), Artem Mamonov (MGH), Jayashree Kalpathy-Cramer (MGH), Andriy Fedorov (BWH)")
set(EXTENSION_DESCRIPTION "Estimate effective T1 from multi-spectral FLASH MRI scans with arbitrary number of flip angles")
set(EXTENSION_ICONURL "http://wiki.slicer.org/slicerWiki/images/3/32/T1_Mapping_Logo_Resized.png")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/4/4c/T1_Mapping_GUI.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(T1_Mapping)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
