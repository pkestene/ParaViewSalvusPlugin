# state file tested using paraview version 5.12
#
# prototype reader for distributed (MPI-based ParaView) reading
# of Salvus meshes.
#
# Written by Jean M. Favre, Swiss National Supercomputing Center
#
# tested on Piz Daint, Wed 13 Dec 10:45:22 CET 2023
#
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

materialLibrary1 = GetMaterialLibrary()

renderView1 = GetRenderView()
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.02780125617980957, -0.027796482086181643, 0.2720464096069336]
renderView1.CameraPosition = [-0.17874834523276456, 0.49986769625918587, 0.18586292891207112]
renderView1.CameraFocalPoint = [-0.02780125617980957, -0.027796482086181643, 0.2720464096069336]
renderView1.CameraViewUp = [0.09418185259355806, 0.18666272160823183, 0.9778991803881693]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.14378838770274488
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# create a new 'Programmable Source'
aCOUSTIC = ProgrammableSource(registrationName='ACOUSTIC')
aCOUSTIC.OutputDataSetType = 'vtkUnstructuredGrid'
aCOUSTIC.Script = """
import h5py
from vtk import vtkUnstructuredGrid, VTK_HEXAHEDRON
from vtkmodules.vtkFiltersCore import vtkStaticCleanUnstructuredGrid
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np

executive = self.GetExecutive()
outInfo = executive.GetOutputInformation(0)
piece = outInfo.Get(executive.UPDATE_PIECE_NUMBER())
numPieces = outInfo.Get(executive.UPDATE_NUMBER_OF_PIECES())

f = h5py.File("/scratch/snx3000/pmarty/tmp/job_2311281715730884_10bd674d21/wavefield.h5","r")

mesh = 'ACOUSTIC'

timestep = 10

dataset1 = 'connectivity_' + mesh
dataset2 = 'coordinates_' + mesh
shape = f[dataset2].shape
total_npts = shape[0]*shape[1]
total_nelts = f[dataset1].shape[0]

if numPieces == 1:
  load = MyNumber_of_Cells = total_nelts
  MyNumber_of_Nodes = total_npts; minId = 0
else:
  load = total_nelts // numPieces
if piece < (numPieces-1):
  MyNumber_of_Cells = load
else:
  MyNumber_of_Cells = total_nelts - (numPieces-1) * load
  
begin = piece*load
end = piece*load + MyNumber_of_Cells

CONNECTIVITY = np.full([MyNumber_of_Cells, 9], 8, dtype='int64')
CONNECTIVITY[:,1:] = f[dataset1][begin:end]
minId = np.min(CONNECTIVITY[:,1:])
maxId = np.max(CONNECTIVITY[:,1:])
if minId != 0:
  CONNECTIVITY[:,1:] = CONNECTIVITY[:,1:] - minId
CELL_TYPES = np.full([MyNumber_of_Cells], VTK_HEXAHEDRON, np.ubyte)

tmpOutput = vtkUnstructuredGrid()
ug = dsa.WrapDataObject(tmpOutput)
ug.SetCells(CELL_TYPES, np.arange(0, MyNumber_of_Cells*9, 9), CONNECTIVITY)

ug.Points = f[dataset2][()].reshape(shape[0]*shape[1], shape[2])[minId:maxId+1:,]
if mesh == 'ACOUSTIC':
  ug.PointData.append(f['volume']['phi_tt'][timestep,:,:,:-3].reshape(shape[0]*shape[1])[minId:maxId+1:,].ravel(), "phi_tt")
else:
  #varnames = ["stress_xx", "stress_yy", "stress_zz", "stress_yz", "stress_xz", "stress_xy"]
  varnames = ["stress_xx", "stress_yy", "stress_zz"]
  for i, var in enumerate(varnames):
    ug.PointData.append(f['volume']['stress'][timestep,:,i,:-3].reshape(shape[0]*shape[1])[minId:maxId+1:,].ravel(), var)
f.close()

cleaner = vtkStaticCleanUnstructuredGrid()
cleaner.SetInputData(tmpOutput)
cleaner.Update()

output.ShallowCopy(cleaner.GetOutput())
"""

# create a new 'Contour'
Acoustic_contours = Contour(registrationName='Acoustic_contours', Input=aCOUSTIC)
Acoustic_contours.ContourBy = ['POINTS', 'phi_tt']
Acoustic_contours.Isosurfaces = [-4e+17, -3.111111111111111e+17, -2.222222222222222e+17, -1.3333333333333331e+17, -4.444444444444442e+16, 4.444444444444448e+16, 1.3333333333333338e+17, 2.2222222222222234e+17, 3.111111111111112e+17, 4e+17]

Acoustic_contoursDisplay = Show(Acoustic_contours)
Acoustic_contoursDisplay.Representation = 'Surface'
ColorBy(Acoustic_contoursDisplay, ['POINTS', 'result'])
