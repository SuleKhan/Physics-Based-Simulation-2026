# GLFW >= 3.4, for Wayland support.
#
# polyscope vendors GLFW 3.3.9, whose Wayland backend is experimental and
# mutually exclusive with X11 at compile time. 3.4 builds both backends and
# chooses at runtime: Wayland when WAYLAND_DISPLAY is set, X11 otherwise, so
# X11-only setups keep working.
#
# polyscope guards its bundled copy with `if(NOT TARGET glfw)`, so defining the
# target before add_subdirectory(ThirdParty) makes polyscope adopt this one.

if(TARGET glfw)
    return()
endif()

include(FetchContent)

set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)

# GLFW_BUILD_X11 and GLFW_BUILD_WAYLAND are cmake_dependent_options that already
# default to ON on Linux and to OFF on macOS and Windows, so both backends are
# built where they exist and neither needs forcing here.

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
)
FetchContent_MakeAvailable(glfw)
