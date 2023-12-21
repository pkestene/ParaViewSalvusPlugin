
set(classes
   vtkSalvusHDF5Reader
   )

set(private_headers
   vtkSalvusHDF5Reader.h
   )
  
vtk_module_add_module(SalvusHDF5Reader
  CLASSES ${classes}
  SOURCES ${sources}
  PRIVATE_HEADERS ${private_headers})
  
paraview_add_server_manager_xmls(XMLS  SalvusHDF5_Server.xml)
