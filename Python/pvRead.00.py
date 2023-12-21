# state file generated using paraview version 5.11.2
# tested on Piz Daint with a single MPI task on a single node
###
#!/bin/bash -l
#SBATCH --job-name=Daint_ParaViewServer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --hint=multithread
#SBATCH --time=00:29:59
#SBATCH --account=csstaff
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#module load daint-gpu
#module unload ddt
#module load h5py
#module load ParaView/5.11.2-CrayGNU-21.09-EGL
#SBATCH --cpus-per-task=16

#srun -n 1 -N 1 -c 16 pvserver --reverse-connection --client-host=daint101.cscs.ch --server-port=110
###
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
LoadPlugin("/users/jfavre/Projects/Salvus/ParaViewSalvusPlugin/build511/./lib64/paraview-5.11/plugins/pvSalvusHDF5Reader/pvSalvusHDF5Reader.so", ns=globals())

reader = SalvusHDF5reader(registrationName='wavefield.h5', FileName='/scratch/snx3000/pmarty/tmp/job_2311281715730884_10bd674d21/wavefield.h5')
reader.PointArrays = ["stress_xx"]
reader.UpdatePipelineInformation()

#reader.GetDataInformation().GetBounds()
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

view = GetRenderView()
view.ViewTime = reader.TimestepValues[-1]

rep = Show(reader)
rep.Representation = 'Outline'
Render()
