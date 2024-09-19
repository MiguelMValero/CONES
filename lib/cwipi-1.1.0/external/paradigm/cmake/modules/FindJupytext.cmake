
find_program(
  JUPYTEXT_EXECUTABLE
  NAMES jupytext
  PATHS ${JUPYTEXT_ROOT} ENV JUPYTEXT_ROOT
  DOC "Path to jupytext executable"
)

if(JUPYTEXT_EXECUTABLE)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    Jupytext DEFAULT_MESSAGE JUPYTEXT_EXECUTABLE
  )
endif()
