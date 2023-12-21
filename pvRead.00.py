# state file generated using paraview version 5.11.2
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
LoadPlugin("/users/jfavre/Projects/Salvus/ParaViewSalvusPlugin/build511/./lib64/paraview-5.11/plugins/pvSalvusHDF5Reader/pvSalvusHDF5Reader.so", ns=globals())

reader = SalvusHDF5reader(registrationName='wavefield.h5', FileName='/scratch/snx3000/pmarty/tmp/job_2311281715730884_10bd674d21/wavefield.h5')
reader.PointArrays = []
reader.UpdatePipelineInformation()

#reader.GetDataInformation().GetBounds()
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

view = GetRenderView()
view.ViewTime = reader.TimestepValues[-1]

rep = Show(reader)
rep.Representation = 'Outline'
Render()
