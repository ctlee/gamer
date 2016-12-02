



set(BLENDER_EXECUTABLE "")
set(BLENDER_VERSION "")
if(APPLE)
    # First look to see if it's installed as an app bundle...
    set(bundle "/Applications/blender.app")
    if(EXISTS ${bundle} AND EXISTS "${bundle}/Contents/Info.plist")
        # Parse Info.plist for exe and version
        set(line_is_main_executable 0)
        set(line_is_version 0)

        set(eol_char "E")
        file(READ "${bundle}/Contents/Info.plist" info_plist)
        string(REPLACE ";" "\\;" info_plist "${info_plist}")
        string(REPLACE "\n" "${eol_char};" info_plist "${info_plist}")
        string(REPLACE "\r" "${eol_char};" info_plist "${info_plist}")

        foreach(line ${info_plist})
            if(line_is_main_executable)
                string(REGEX REPLACE "^.*<string>(.*)</string>.*$" "\\1" bundle_executable "${line}")
                set(line_is_main_executable 0)
            elseif(line_is_version)
                string(REGEX REPLACE "^.*<string>(.*)</string>.*$" "\\1" BLENDER_VERSION "${line}")
                break()
            endif()

            if(line MATCHES "<key>CFBundleExecutable</key>")
                set(line_is_main_executable 1)
            elseif(line MATCHES "<key>CFBundleShortVersionString</key>")
                set(line_is_version 1)
            endif()
        endforeach()
        if(NOT "${bundle_executable}" STREQUAL "")
            if(EXISTS "${bundle}/Contents/MacOS/${bundle_executable}")
                set(BLENDER_EXECUTABLE "${bundle}/Contents/MacOS/${bundle_executable}")
            endif()
        endif()
        if(NOT "${BLENDER_VERSION}" STREQUAL "")
            if(EXISTS "${bundle}/Contents/Resources/${BLENDER_VERSION}/python")
            endif()
        endif()
    endif()
elseif(WIN32)
    # TODO figure out how to best search for windows bin
endif()

# If we haven't found the exe, do a cmake search
if("${BLENDER_EXECUTABLE}" STREQUAL "")
    find_program(blender_exe
            blender)
    if(blender_exe)
        set(BLENDER_EXECUTABLE ${blender_exe})
        # Parse for the version by running blender
        execute_process(COMMAND ${BLENDER_EXECUTABLE} -v
                        OUTPUT_VARIABLE blender_version_long
                        ERROR_QUIET
                        OUTPUT_STRIP_TRAILING_WHITESPACE) 
        string(REGEX REPLACE "^.*Blender ([0-9]*.[0-9]*) \(.*\).*$" "\\1" BLENDER_VERSION "${blender_version_long}")
    else()
        message(FATAL_ERROR "Blender not found!")
    endif()
endif()

message(STATUS "Blender executable: ${BLENDER_EXECUTABLE}")
message(STATUS "Blender version: ${BLENDER_VERSION}")

if(BLENDER_EXECUTABLE)
    # Get the python path and version from blender
    # This is kind of a hack... Is there a better way?
    execute_process(COMMAND ${BLENDER_EXECUTABLE} -b --factory-startup --python-expr "import sys;print(sys.path[1]);print('zz' + ';'.join([str(x) for x in sys.version_info[:3]]) + 'zz')"
                    OUTPUT_VARIABLE blender_output
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    
    set(eol_char "E")
    string(REPLACE "\n" "${eol_char};" blender_output "${blender_output}")
    string(REPLACE "\r" "${eol_char};" blender_output "${blender_output}")
    string(REGEX REPLACE "^.*found bundled python: (.*python).*$" "\\1" BLENDER_PYTHON_ROOT "${blender_output}")
    string(REGEX REPLACE "(.*)/addons.*$" "\\1" BLENDER_ADDON_PATH "${blender_output}")
    string(REGEX REPLACE "^.*zz(.*)zz.*$" "\\1" _VERSION "${blender_output}")
    
    if(_VERSION)
        string(REPLACE ";" "." BLENDER_PYTHON_VERSION_STRING "${_VERSION}")
        list(GET _VERSION 0 BLENDER_PYTHON_VERSION_MAJOR)
        list(GET _VERSION 1 BLENDER_PYTHON_VERSION_MINOR)
        list(GET _VERSION 2 BLENDER_PYTHON_VERSION_PATCH)
        if(PYTHON_VERSION_PATCH EQUAL 0)
            # it's called "Python 2.7", not "2.7.0"
            string(REGEX REPLACE "\\.0$" "" BLENDER_PYTHON_VERSION_STRING "${BLENDER_PYTHON_VERSION_STRING}")
        endif()
    endif()
endif()

message(STATUS "Blender python path: ${BLENDER_PYTHON_ROOT}")
message(STATUS "Blender python version: ${BLENDER_PYTHON_VERSION_STRING}")
message(STATUS "Blender addon path: ${BLENDER_ADDON_PATH}")
