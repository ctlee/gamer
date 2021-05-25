
find_package(Python QUIET COMPONENTS Interpreter )
if(Python_Interpreter_FOUND)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "from importlib import util;print(util.find_spec('exhale') is not None)"
                    OUTPUT_VARIABLE EXHALE_IMPORT_STATUS
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(Python_Interpreter_FOUND)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Exhale DEFAULT_MSG
    EXHALE_IMPORT_STATUS
)
