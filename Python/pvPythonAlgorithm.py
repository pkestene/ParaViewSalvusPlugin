# plugin tested using paraview version 5.12
#
# prototype reader for distributed (MPI-based ParaView) reading
# of Salvus meshes.
#
# Written by Jean M. Favre, Swiss National Supercomputing Center
#
# tested on Piz Daint, Thu 11 Jan 06:31:10 CET 2024
#
from paraview.util.vtkAlgorithm import *
import numpy as np
import h5py
from vtkmodules.numpy_interface import dataset_adapter as dsa

def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

@smproxy.reader(name="SalvusHDF5Reader", label="Python Salvus HDF5 Reader", extensions="h5", file_description="HDF5 files")
class SalvusHDF5Reader(VTKPythonAlgorithmBase):
    """A reader that reads a Salvus HDF5 file."""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkUnstructuredGrid')
        self._filename = None
        self._model = None
        self._ndata = None
        self._timesteps = None
        self.dt = 1.0
        self.t_start = 0.0
        self.cached_geometry = None
        self.total_npts = -1
        self.total_nelts = -1
        self.minId = -1
        self.maxId = -1

        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))

    def _get_timesteps(self):
        if self._filename is not None:
          if self._timesteps is None:
            with h5py.File(self._filename) as f:
              sampling_rate = f["volume"].attrs['sampling_rate_in_hertz']
              self.t_start = f["volume"].attrs['start_time_in_seconds']
              if 'phi_tt' in f["volume"].keys():
                NumLoadSteps = f["volume"]["phi_tt"].shape[0]
              elif 'stress' in f["volume"].keys():
                NumLoadSteps = f["volume"]["stress"].shape[0]
              self.dt = 1.0 /  sampling_rate
              self._timesteps = np.arange(NumLoadSteps) * self.dt + self.t_start
              print("Read timesteps from file: ", self._timesteps)
        return self._timesteps.tolist() if self._timesteps is not None else None

    def _get_update_piece(self, outInfo):
        executive = self.GetExecutive()
        piece = outInfo.Get(executive.UPDATE_PIECE_NUMBER())
        numPieces = outInfo.Get(executive.UPDATE_NUMBER_OF_PIECES())
        
        if numPieces == 1:
          load = MyNumber_of_Cells = self.total_nelts
          MyNumber_of_Nodes = self.total_nelts; minId = 0
        else:
          load = self.total_nelts // numPieces
        if piece < (numPieces-1):
          MyNumber_of_Cells = load
        else:
          MyNumber_of_Cells = self.total_nelts - (numPieces-1) * load
  
        begin = piece*load
        end = piece*load + MyNumber_of_Cells
        return MyNumber_of_Cells, begin, end
            
    def _get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self._get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    def _get_array_selection(self):
        return self._arrayselection

    # "StringInfo" and "String" demonstrate how one can add a selection widget
    # that lets user choose a string from the list of strings.
    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return ["ELASTIC", "ACOUSTIC"]

    @smproperty.stringvector(name="Model", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="StringInfo" function="StringInfo"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetModel(self, value):
        """Specify model to file (ACOUSTIC or ELASTIC"""
        if value in self.GetStrings() and self._model != value:
            print("\nModel ", value, " was selected\n")
            self._model = value
            self.Modified()

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="h5", file_description="Salvus HDF5 files")
    def SetFileName(self, name):
        if self._filename != name:
            self._filename = name
            self._ndata = None
            self._timesteps = None
            self.Modified()

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._get_timesteps()

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonExecutionModel import vtkAlgorithm
        print("\nRequestInformation: begin\n")
        
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Set(vtkAlgorithm.CAN_HANDLE_PIECE_REQUEST(), 1)

        if self._model in self.GetStrings():
            with h5py.File(self._filename) as f:
                dataset1 = 'connectivity_' + self._model
                dataset2 = 'coordinates_' + self._model
                shape = f[dataset2].shape
                self.total_npts = shape[0]*shape[1]
                self.total_nelts = f[dataset1].shape[0]
            
        for aname in ["phi_tt","stress_xx", "stress_yy", "stress_zz", "stress_yz", "stress_xz", "stress_xy"]:
            self._arrayselection.AddArray(aname)
            
        timesteps = self._get_timesteps()
        if timesteps is not None:
            outInfo.Remove(executive.TIME_STEPS())
            outInfo.Remove(executive.TIME_RANGE())
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])

        print("\n===============================================: end\n")
        return 1

    def Create_new_mesh(self, outInfoVec, t_index):
        from vtk import VTK_HEXAHEDRON, vtkUnstructuredGrid
        MyNumber_of_Cells, begin, end = self._get_update_piece(outInfoVec.GetInformationObject(0))
        #create initial instance of geometry
        ug0 = vtkUnstructuredGrid()
        out0 = dsa.WrapDataObject(ug0)
        
        with h5py.File(self._filename) as f:
          CONNECTIVITY = np.full([MyNumber_of_Cells, 9], 8, dtype='int64')
          CONNECTIVITY[:,1:] = f['connectivity_' + self._model][begin:end]
          self.minId = np.min(CONNECTIVITY[:,1:])
          self.maxId = np.max(CONNECTIVITY[:,1:])
          dataset2 = 'coordinates_' + self._model
          shape = f[dataset2].shape
          if self.minId != 0:
            CONNECTIVITY[:,1:] = CONNECTIVITY[:,1:] - self.minId
          CELL_TYPES = np.full([MyNumber_of_Cells], VTK_HEXAHEDRON, np.ubyte)
          out0.SetCells(CELL_TYPES, np.arange(0, MyNumber_of_Cells*9, 9), CONNECTIVITY)
          out0.Points = f[dataset2][()].reshape(shape[0]*shape[1], shape[2])[self.minId:self.maxId+1:,]

          narrays = self._arrayselection.GetNumberOfArrays()
          for i in range(narrays):
            if self._arrayselection.GetArraySetting(i):
              name = self._arrayselection.GetArrayName(i) # using (i-1) for ELASTIC. TODO. FIXME for ACOUSTIC
              out0.PointData.append(f['volume']['stress'][t_index,:,i-1,:-3].reshape(shape[0]*shape[1])[self.minId:self.maxId+1:,].ravel(), name)

        return ug0
        
    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        from vtk import VTK_HEXAHEDRON

        data_time = self._get_update_time(outInfoVec.GetInformationObject(0))
        # need now to find my index in array timesteps
        t_index = np.where(self._timesteps == data_time)[0][0]
        print("RequestData: self._timesteps[", t_index, "] = ", data_time)

        ug = vtkUnstructuredGrid.GetData(outInfoVec, 0)
        output = dsa.WrapDataObject(ug)
        if self.cached_geometry is None:
          ug1 = self.Create_new_mesh(outInfoVec, t_index)
          print("creating a new geometry & saving it in cache")
          self.cached_geometry = ug1
          ug.ShallowCopy(ug1)
        else:
          ug.ShallowCopy(self.cached_geometry)
          print("reusing the cached geometry")
          narrays = self._arrayselection.GetNumberOfArrays()
          with h5py.File(self._filename) as f:
            dataset2 = 'coordinates_' + self._model
            shape = f[dataset2].shape
            for i in range(narrays):
              if self._arrayselection.GetArraySetting(i):
                name = self._arrayselection.GetArrayName(i) # using (i-1) for ELASTIC. TODO. FIXME for ACOUSTIC
                output.PointData.append(f['volume']['stress'][t_index,:,i-1,:-3].reshape(shape[0]*shape[1])[self.minId:self.maxId+1:,].ravel(), name)

        if data_time is not None:
            output.GetInformation().Set(output.DATA_TIME_STEP(), data_time)
        print("\n==================: end\n")
        return 1

def test_SalvusHDF5Reader(fname):
    reader = SalvusHDF5Reader()
    reader.SetFileName(fname)
    reader.Update()
    assert reader.GetOutputDataObject(0).GetNumberOfPoints() > 0

if __name__ == "__main__":
    print("testing...testing...testing")
    test_SalvusHDF5Reader("/local/data/Salvus/wavefield.h5")

    from paraview.detail.pythonalgorithm import get_plugin_xmls
    from xml.dom.minidom import parseString
    for xml in get_plugin_xmls(globals()):
        dom = parseString(xml)
        print(dom.toprettyxml(" ","\n"))

