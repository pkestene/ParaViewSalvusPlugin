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
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.02779959297180176, -0.027796482086181643, 0.2720464096069336]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.17831580106353462, 0.3491437129910117, 0.39797301762884874]
renderView1.CameraFocalPoint = [0.016051944274902447, -0.1418534561528596, 0.2253996955506869]
renderView1.CameraViewUp = [0.1919441785852059, -0.2569259728701931, 0.9471781652740712]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.1437875702964892
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# create a new 'XML PolyData Reader'
skullvtp = XMLPolyDataReader(registrationName='skull.vtp', FileName=['/users/jfavre/skull.vtp'])
skullvtp.PointArrayStatus = ['stress_xx', 'stress_yy', 'stress_zz']
skullvtp.TimeArray = 'None'

# create a new 'Python Calculator'
pythonCalculator2 = PythonCalculator(registrationName='PythonCalculator2', Input=skullvtp)
pythonCalculator2.Expression = '(stress_xx+stress_yy+stress_zz)/3.0'

# create a new 'Programmable Source'
eLASTIC = ProgrammableSource(registrationName='ELASTIC')
eLASTIC.OutputDataSetType = 'vtkUnstructuredGrid'
eLASTIC.Script = """
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

mesh = 'ELASTIC'

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

CONNECTIVITY = np.full([MyNumber_of_Cells, 9], 8, dtype=\'int64\')
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
  ug.PointData.append(f['volume\']['phi_tt'][timestep,:,:,:-3].reshape(shape[0]*shape[1])[minId:maxId+1:,].ravel(), "phi_tt")
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

# create a new 'Clip'
Elastic_clip1 = Clip(registrationName='Elastic_clip1', Input=pythonCalculator2)
Elastic_clip1.ClipType = 'Plane'
Elastic_clip1.Invert = 0
Elastic_clip1.ClipType.Origin = [-0.02780125617980957, -0.027796482086181643, 0.2720464096069336]

# create a new 'Clip'
Elastic_clip3 = Clip(registrationName='Elastic_clip3', Input=Elastic_clip1)
Elastic_clip3.ClipType = 'Plane'
Elastic_clip3.Invert = 0
Elastic_clip3.ClipType.Origin = [0.007532441139221191, -0.02775630950927735, 0.2722850952148438]
Elastic_clip3.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Clip'
Elastic_clip2 = Clip(registrationName='Elastic_clip2', Input=pythonCalculator2)
Elastic_clip2.ClipType = 'Plane'
Elastic_clip2.ClipType.Origin = [-0.02780125617980957, -0.027796482086181643, 0.2720464096069336]
Elastic_clip2.ClipType.Normal = [0.0, 0.0, 1.0]

#Elastic_cleantoGrid = CleantoGrid(registrationName='Elastic_cleantoGrid', Input=eLASTIC)

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(registrationName='PythonCalculator1', Input=eLASTIC)
pythonCalculator1.Expression = '(stress_xx+stress_yy+stress_zz)/3.0'

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=skullvtp)

# create a new 'Contour'
Elastic_contours = Contour(registrationName='Elastic_contours', Input=pythonCalculator1)
Elastic_contours.ContourBy = ['POINTS', 'result']
Elastic_contours.Isosurfaces = [-4e+17, -3.111111111111111e+17, -2.222222222222222e+17, -1.3333333333333331e+17, -4.444444444444442e+16, 4.444444444444448e+16, 1.3333333333333338e+17, 2.2222222222222234e+17, 3.111111111111112e+17, 4e+17]

Elastic_contoursDisplay = Show(Elastic_contours)

# get color transfer function/color map for 'result'
resultLUT = GetColorTransferFunction('result')
resultLUT.RGBPoints = [-4.697729056773243e+17, 0.231373, 0.298039, 0.752941, -2.693623099424768e+16, 0.865003, 0.865003, 0.865003, 4.159004436888289e+17, 0.705882, 0.0156863, 0.14902]
resultLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
Elastic_contoursDisplay.Representation = 'Surface'
Elastic_contoursDisplay.ColorArrayName = ['POINTS', 'result']
Elastic_contoursDisplay.LookupTable = resultLUT

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
Elastic_contoursDisplay.ScaleTransferFunction.Points = [-4e+17, 0.0, 0.5, 0.0, 4e+17, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
Elastic_contoursDisplay.OpacityTransferFunction.Points = [-4e+17, 0.0, 0.5, 0.0, 4e+17, 1.0, 0.5, 0.0]

# show data from Elastic_clip2
Elastic_clip2Display = Show(Elastic_clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
Elastic_clip2Display.Representation = 'Surface'
Elastic_clip2Display.ColorArrayName = ['POINTS', '']
Elastic_clip2Display.SelectTCoordArray = 'None'
Elastic_clip2Display.SelectNormalArray = 'Normals'
Elastic_clip2Display.SelectTangentArray = 'None'
Elastic_clip2Display.BackfaceRepresentation = 'Surface'
Elastic_clip2Display.BackfaceDiffuseColor = [0.666666, 0.666666, 0.666666]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
Elastic_clip2Display.ScaleTransferFunction.Points = [-0.9999999403953552, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
Elastic_clip2Display.OpacityTransferFunction.Points = [-0.9999999403953552, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from Elastic_clip3
Elastic_clip3Display = Show(Elastic_clip3, renderView1, 'UnstructuredGridRepresentation')

# get opacity transfer function/opacity map for 'result'
resultPWF = GetOpacityTransferFunction('result')
resultPWF.Points = [-4.697729056773243e+17, 0.0, 0.5, 0.0, 4.159004436888289e+17, 1.0, 0.5, 0.0]
resultPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
Elastic_clip3Display.Representation = 'Surface'
Elastic_clip3Display.ColorArrayName = ['POINTS', 'result']
Elastic_clip3Display.LookupTable = resultLUT
Elastic_clip3Display.SelectTCoordArray = 'None'
Elastic_clip3Display.SelectNormalArray = 'Normals'
Elastic_clip3Display.SelectTangentArray = 'None'
Elastic_clip3Display.BackfaceRepresentation = 'Surface'
Elastic_clip3Display.BackfaceDiffuseColor = [0.666666, 0.666666, 0.666666]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
Elastic_clip3Display.ScaleTransferFunction.Points = [-0.9999999403953552, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
Elastic_clip3Display.OpacityTransferFunction.Points = [-0.9999999403953552, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# get color legend/bar for resultLUT in view renderView1
resultLUTColorBar = GetScalarBar(resultLUT, renderView1)
resultLUTColorBar.Title = 'result'
resultLUTColorBar.ComponentTitle = ''

# set color bar visibility
resultLUTColorBar.Visibility = 1

# show color legend
Elastic_contoursDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
Elastic_clip3Display.SetScalarBarVisibility(renderView1, True)

SetActiveSource(eLASTIC)

