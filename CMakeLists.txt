cmake_minimum_required(VERSION 3.20...4.0)

project(PBS
    DESCRIPTION "Physically-based Simulation in Computer Graphics - course exercises"
    LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Default to an optimized build when no build type is given. Multi-config
# generators (Visual Studio, Xcode) choose per build instead, so leave those
# alone.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

# libigl 2.6.4 still calls FetchContent_Populate() internally to fetch Eigen,
# which CMP0169 (CMake 3.30) promotes from a warning to a hard error. Opt out
# for libigl's scope only; remove once libigl moves to FetchContent_MakeAvailable.
cmake_policy(PUSH)
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()
include(libigl)                 # FetchContent -> igl::core
cmake_policy(POP)
include(glfw)                   # GLFW 3.4 -> glfw, before polyscope can vendor 3.3.9
add_subdirectory(ThirdParty)    # -> polyscope (brings its own glad/imgui)

# Optional on purpose: AppleClang ships no OpenMP runtime, so macOS degrades
# to a serial build rather than failing to configure.
find_package(OpenMP)

# ---------------------------------------------------------------------------
# pbs_common - framework code shared by every exercise
# ---------------------------------------------------------------------------
add_library(pbs_common STATIC
    include/CameraHelper.cpp
    include/GeometryHelper.cpp
    include/Optimization.cpp
    include/PBSApp.cpp)

target_include_directories(pbs_common PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# C++20 for <numbers>. Set on the target rather than via CMAKE_CXX_STANDARD,
# which would also rebuild polyscope, glfw, imgui and glm as C++20.
target_compile_features(pbs_common PUBLIC cxx_std_20)

target_compile_definitions(pbs_common PUBLIC IGL_VIEWER_VIEWER_QUIET)

target_link_libraries(pbs_common PUBLIC igl::core polyscope)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pbs_common PUBLIC OpenMP::OpenMP_CXX)
endif()

# PUBLIC so the exercise targets inherit it. MSVC has no -Wall, and its /Wall
# is unusably noisy.
target_compile_options(pbs_common PUBLIC
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>
    $<$<CXX_COMPILER_ID:MSVC>:/W3>)

# ---------------------------------------------------------------------------
include(PBSExercise)

add_subdirectory(1_mass_spring)
