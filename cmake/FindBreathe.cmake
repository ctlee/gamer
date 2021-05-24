
find_package(Python QUIET COMPONENTS Interpreter )
if(Python_Interpreter_FOUND)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "from importlib import util;print(util.find_spec('breathe') is not None)"
                    OUTPUT_VARIABLE BREATHE_IMPORT_STATUS
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(Python_Interpreter_FOUND)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Breathe DEFAULT_MSG
    BREATHE_IMPORT_STATUS
)
