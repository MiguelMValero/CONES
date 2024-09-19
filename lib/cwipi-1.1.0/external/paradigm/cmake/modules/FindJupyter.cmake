
find_program(
  JUPYTER_EXECUTABLE
  NAMES jupyter
  PATHS ${JUPYTER_ROOT} ENV JUPYTER_ROOT
  DOC "Path to jupyter executable"
)

if(JUPYTER_EXECUTABLE)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    Jupyter DEFAULT_MESSAGE JUPYTER_EXECUTABLE
  )
endif()
