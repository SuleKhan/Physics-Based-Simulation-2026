# add_pbs_exercise(<name> [EXTRA_DIRS <dir>...])
#
# Builds the two executables for one exercise from a shared source layout:
#
#   <name>/*.cpp            compiled into both variants (main.cpp, Test.cpp, ...)
#   <name>/template/*.cpp   the version to fill in    -> target <name>
#   <name>/solution/*.cpp   the reference version     -> target <name>_solution
#
# EXTRA_DIRS names further source directories compiled into both variants;
# 3_fluid keeps its grid helpers in helper/.

function(add_pbs_exercise name)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "" "EXTRA_DIRS")

    foreach(variant IN ITEMS template solution)
        # The handed-out repository ships template/ without solution/, so build
        # only the variants whose directory is actually present.
        if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${variant}")
            continue()
        endif()

        if(variant STREQUAL "template")
            set(target "${name}")
        else()
            set(target "${name}_${variant}")
        endif()

        set(dirs "" "${variant}" ${ARG_EXTRA_DIRS})
        set(patterns "")
        foreach(dir IN LISTS dirs)
            if(dir STREQUAL "")
                list(APPEND patterns "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
                                     "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
            else()
                list(APPEND patterns "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp"
                                     "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.h")
            endif()
        endforeach()

        # CONFIGURE_DEPENDS re-globs at build time, so a newly added source
        # file is picked up without re-running cmake by hand.
        file(GLOB sources CONFIGURE_DEPENDS ${patterns})

        # Headers are listed as sources purely so Visual Studio and Xcode show
        # them in the project tree; they do not affect the compile.
        add_executable(${target} ${sources})
        target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
        target_link_libraries(${target} PRIVATE pbs_common)

        # Absolute path to this exercise's data, fixed at configure time. C++
        # has no portable way to locate the running executable, and a relative
        # path cannot work across both single- and multi-config generators:
        # Visual Studio and Xcode put binaries in a per-config subdirectory and
        # default the debugger's working directory elsewhere again.
        target_compile_definitions(${target}
            PRIVATE PBS_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    endforeach()
endfunction()
