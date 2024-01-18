# state file generated using paraview version 5.12.0-RC2
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [640, 720]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.027801258489489555, -0.027796484529972076, 0.27204640209674835]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.3764026201618533, 0.33594493172346757, 0.5061646903381758]
renderView1.CameraFocalPoint = [-0.027801258489489472, -0.02779648452997201, 0.2720464020967484]
renderView1.CameraViewUp = [0.34758199981765503, -0.24876609976277442, 0.9040476652320812]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.14378838746977188
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [640, 720]
renderView2.AxesGrid = 'Grid Axes 3D Actor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [-0.027801258489489555, -0.027796484529972076, 0.27204640209674835]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-0.3764026201618533, 0.33594493172346757, 0.5061646903381758]
renderView2.CameraFocalPoint = [-0.027801258489489472, -0.02779648452997201, 0.2720464020967484]
renderView2.CameraViewUp = [0.34758199981765503, -0.24876609976277442, 0.9040476652320812]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 0.14378838746977188
renderView2.LegendGrid = 'Legend Grid Actor'
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)
layout1.SetSize(1303, 784)

SetActiveView(renderView1)

fname = '/scratch/snx3000/pmarty/tmp/job_2311281715730884_10bd674d21/wavefield.h5'

# create a new 'Salvus HDF5 reader'
acoustic = SalvusHDF5reader(registrationName='wavefield.h5', FileName=fname)
acoustic.ElasticPointArrays = []
acoustic.AcousticPointArrays = ['phi_tt']
acoustic.ModelName = 'ACOUSTIC'

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=acoustic)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-0.028819311410188675, -0.021390166133642197, 0.282152459025383]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [-0.028819311410188675, -0.021390166133642197, 0.282152459025383]

# create a new 'Salvus HDF5 reader'
elastic = SalvusHDF5reader(registrationName='wavefield.h5', FileName=fname)
elastic.ElasticPointArrays = ['stress_yy']
elastic.AcousticPointArrays = ['phi_tt']

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=elastic)
annotateTimeFilter1.Format = 'Time: {time:g}'

wavefieldh5_1Display = Show(elastic, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'stress_yy'
stress_yyLUT = GetColorTransferFunction('stress_yy')
stress_yyLUT.AutomaticRescaleRangeMode = 'Never'
stress_yyLUT.RGBPoints = [-3e+17, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 3e+17, 0.705882, 0.0156863, 0.14902]
stress_yyLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'stress_yy'
stress_yyPWF = GetOpacityTransferFunction('stress_yy')
stress_yyPWF.Points = [-3e+17, 0.0, 0.5, 0.0, 3e+17, 1.0, 0.5, 0.0]
stress_yyPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
wavefieldh5_1Display.Representation = 'Surface'
wavefieldh5_1Display.ColorArrayName = ['POINTS', 'stress_yy']
wavefieldh5_1Display.LookupTable = stress_yyLUT

wavefieldh5_1Display.OSPRayScaleArray = 'GlobalNodeIds'
wavefieldh5_1Display.OSPRayScaleFunction = 'Piecewise Function'
wavefieldh5_1Display.Assembly = ''
wavefieldh5_1Display.SelectOrientationVectors = 'None'
wavefieldh5_1Display.ScaleFactor = 0.01891310811042786
wavefieldh5_1Display.SelectScaleArray = 'GlobalNodeIds'
wavefieldh5_1Display.GlyphType = 'Arrow'
wavefieldh5_1Display.GlyphTableIndexArray = 'GlobalNodeIds'
wavefieldh5_1Display.GaussianRadius = 0.0009456554055213928
wavefieldh5_1Display.SetScaleArray = ['POINTS', 'GlobalNodeIds']
wavefieldh5_1Display.ScaleTransferFunction = 'Piecewise Function'
wavefieldh5_1Display.OpacityArray = ['POINTS', 'GlobalNodeIds']
wavefieldh5_1Display.OpacityTransferFunction = 'Piecewise Function'
wavefieldh5_1Display.DataAxesGrid = 'Grid Axes Representation'
wavefieldh5_1Display.PolarAxes = 'Polar Axes Representation'
wavefieldh5_1Display.ScalarOpacityFunction = stress_yyPWF
wavefieldh5_1Display.ScalarOpacityUnitDistance = 0.0006967230779772128
wavefieldh5_1Display.OpacityArrayName = ['POINTS', 'GlobalNodeIds']

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
wavefieldh5_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 137344374.0, 1.0, 0.5, 0.0]
wavefieldh5_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 137344374.0, 1.0, 0.5, 0.0]

# show data from annotateTimeFilter1
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# get color legend/bar for stress_yyLUT in view renderView1
stress_yyLUTColorBar = GetScalarBar(stress_yyLUT, renderView1)
stress_yyLUTColorBar.Orientation = 'Horizontal'
stress_yyLUTColorBar.WindowLocation = 'Any Location'
stress_yyLUTColorBar.Position = [0.55, 0.08]
stress_yyLUTColorBar.Title = 'stress_yy'
stress_yyLUTColorBar.ComponentTitle = ''
stress_yyLUTColorBar.ScalarBarLength = 0.33

# set color bar visibility
stress_yyLUTColorBar.Visibility = 1

# show color legend
wavefieldh5_1Display.SetScalarBarVisibility(renderView1, True)

slice1Display = Show(slice1, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'phi_tt'
phi_ttLUT = GetColorTransferFunction('phi_tt')
phi_ttLUT.AutomaticRescaleRangeMode = 'Never'

phi_ttLUT.RGBPoints = [-1e+18, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 1e+18, 0.705882, 0.0156863, 0.14902]
phi_ttLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'phi_tt']
slice1Display.LookupTable = phi_ttLUT

slice1Display.OSPRayScaleArray = 'phi_tt'
slice1Display.OSPRayScaleFunction = 'Piecewise Function'
slice1Display.Assembly = ''
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.022280994802713394
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.0011140497401356698
slice1Display.SetScaleArray = ['POINTS', 'phi_tt']
slice1Display.ScaleTransferFunction = 'Piecewise Function'
slice1Display.OpacityArray = ['POINTS', 'phi_tt']
slice1Display.OpacityTransferFunction = 'Piecewise Function'
slice1Display.DataAxesGrid = 'Grid Axes Representation'
slice1Display.PolarAxes = 'Polar Axes Representation'


# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-3.5013868071933706e+18, 0.0, 0.5, 0.0, 2.3007727487811584e+18, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-3.5013868071933706e+18, 0.0, 0.5, 0.0, 2.3007727487811584e+18, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for phi_ttLUT in view renderView2
phi_ttLUTColorBar = GetScalarBar(phi_ttLUT, renderView2)
phi_ttLUTColorBar.Orientation = 'Horizontal'
phi_ttLUTColorBar.WindowLocation = 'Any Location'
phi_ttLUTColorBar.Position = [0.55, 0.08]
phi_ttLUTColorBar.Title = 'phi_tt'
phi_ttLUTColorBar.ComponentTitle = ''
phi_ttLUTColorBar.ScalarBarLength = 0.33

# set color bar visibility
phi_ttLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView2, True)

# get opacity transfer function/opacity map for 'phi_tt'
phi_ttPWF = GetOpacityTransferFunction('phi_tt')
phi_ttPWF.Points = [-1e+18, 0.0, 0.5, 0.0, 1e+18, 1.0, 0.5, 0.0]
phi_ttPWF.ScalarRangeInitialized = 1

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = [renderView1, renderView2]
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 4.289801719272077e-05
animationScene1.StartTime = -6.237574409869409e-06
animationScene1.EndTime = 4.289801719272077e-05
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(elastic)
# ----------------------------------------------------------------

SaveAnimation('/users/jfavre/Projects/Salvus/animation.png', layout1, 16, SaveAllViews=1,
    ImageResolution=[1280, 720],
    FrameWindow=[0, 10])

