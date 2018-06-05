

set(NUMPY_FOUND "")

# Take hints from pythonlibs and pythoninterp
if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND ${PYTHON_VERSION_STRING} STREQUAL ${PYTHONLIBS_VERSION_STRING})
    string(REGEX REPLACE "^(/.*)/.*$" "\\1" PYTHON_LIBS_ROOT ${PYTHON_LIBRARIES})
    message(STATUS "Python site-packages root: ${PYTHON_LIBS_ROOT}")

set(python_sp_root "${PYTHON_LIBS_ROOT}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
if(NOT EXISTS "${python_sp_root}")
        message(SEND_ERROR "Python site-packages could not be found.") 
    else()
        set(numpy_root "${python_sp_root}/numpy")
        if(EXISTS "${numpy_root}")
            set(NUMPY_FOUND "1")
            if(EXISTS "${numpy_root}/core/include")
                set(NUMPY_INCLUDE_DIRS "${numpy_root}/core/include")
            else()
                message(SEND_ERROR "Numpy include dir couldn't be found")
            endif()

            # Read verion.py for short-version
            if(EXISTS "${numpy_root}/version.py")
                file(READ "${numpy_root}/version.py" np_version_file)
                set(eol_char "E")
                string(REPLACE ";" "\\;" np_version_file "${np_version_file}")
                string(REPLACE "\n" "${eol_char};" np_version_file "${np_version_file}")
                string(REPLACE "\r" "${eol_char};" np_version_file "${np_version_file}")
                
                foreach(line ${np_version_file})  
                    if(line MATCHES "short_version.*'")
                        string(REGEX REPLACE "^.*short_version = '(.*)'.*$" "\\1" NUMPY_VERSION ${line})
                        #message(STATUS "Numpy version is: ${NUMPY_VERSION}")
                        break()
                    endif()
                endforeach()
                # TODO: Split the np_version string into components
            elseif()
                # Some other way of getting np version?
            endif()
        else()
            # Numpy doesn't exist here...
        endif()
    endif()
endif()

