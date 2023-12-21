paraview_add_plugin(pvSalvusHDF5Reader
  VERSION "1.0"
  REQUIRED_ON_SERVER
  MODULES SalvusHDF5Reader
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")

option(BUILD_TESTING "Build Plugin Testing" OFF)
if (BUILD_TESTING AND BUILD_SHARED_LIBS)
  add_subdirectory(Testing)
endif()
