
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from importlib import util;print(util.find_spec('breathe') is not None)"
                    OUTPUT_VARIABLE BREATHE_IMPORT_STATUS
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(PYTHONINTERP_FOUND)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(breathe DEFAULT_MSG
    BREATHE_IMPORT_STATUS
)