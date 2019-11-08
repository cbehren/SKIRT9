# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SKIRT project
# ------------------------------------------------------------------

# define project and ensure minimum level of functionality
project(SKIRTproject)
cmake_minimum_required(VERSION 3.2.2 FATAL_ERROR)

# define a user-configurable option to build SKIRT
option(BUILD_SKIRT "build SKIRT, advanced radiative transfer" ON)

# define a user-configurable option to build MakeUp, which requires Qt5
option(BUILD_MAKE_UP "build MakeUp, desktop GUI wizard - requires Qt5")

# define a user-configurable option to build doxstyle
option(BUILD_DOX_STYLE "build doxstyle, documentation block streamliner")

# optionally require HDF5
option(ENABLE_HDF5 "Enable HDF5 support" ON)
if(ENABLE_HDF5)
    include_directories(../../HighFive)
    find_package(HDF5)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs hdf5 hdf5_cpp)
        message("Enabling HDF5 HighFive library support")
        find_package(HighFive 2.0 REQUIRED)
        include_directories($<TARGET_PROPERTY:HighFive,INTERFACE_INCLUDE_DIRECTORIES>)
        link_libraries($<TARGET_PROPERTY:HighFive,INTERFACE_LINK_LIBRARIES>)
    endif()
endif()
# add all relevant subdirectories; each subdirectory defines a single target
add_subdirectory(SMILE)
if (BUILD_SKIRT)
    add_subdirectory(SKIRT)
endif()
if (BUILD_MAKE_UP)
    add_subdirectory(MakeUp)
endif()
if (BUILD_DOX_STYLE)
    add_subdirectory(Docs/doxstyle)
endif()

# set the build type to "Release" as a default; this can be overridden by the user
set(CMAKE_BUILD_TYPE Release)

# mark some obscure variables as advanced so they don't clutter the user interface
mark_as_advanced(CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE CMAKE_CODEBLOCKS_MAKE_ARGUMENTS QT_QMAKE_EXECUTABLE)
mark_as_advanced(CORE_FOUNDATION)
