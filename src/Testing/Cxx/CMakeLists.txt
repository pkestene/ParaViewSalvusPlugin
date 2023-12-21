
ADD_EXECUTABLE(TestSalvusHDF5Reader TestSalvusHDF5Reader.cxx)
TARGET_LINK_LIBRARIES(TestSalvusHDF5Reader
	PUBLIC SalvusHDF5Reader
	PRIVATE
	  VTK::CommonCore
	  VTK::CommonDataModel
          )

