
find_program(
  BREATHE_APIDOC_EXECUTABLE
  NAMES breathe-apidoc
  PATHS ${BREATHE_APIDOC_ROOT} ENV BREATHE_APIDOC_ROOT
  DOC "Path to breathe-apidoc executable"
)

if(BREATHE_APIDOC_EXECUTABLE)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    Breathe DEFAULT_MESSAGE BREATHE_APIDOC_EXECUTABLE
  )
endif()
