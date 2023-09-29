# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
cycle_1txt = CSVReader(registrationName='cycle_1.txt*', FileName=['/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/output/cycle_1.txt', '/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/output/cycle_2.txt'])

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
cycle_1txtDisplay = Show(cycle_1txt, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# Properties modified on cycle_1txtDisplay
cycle_1txtDisplay.Assembly = ''

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=cycle_1txt)
tableToPoints1.XColumn = 'TSL_norm'
tableToPoints1.YColumn = 'TSL_norm'
tableToPoints1.ZColumn = 'TSL_norm'

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'y'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.SelectTCoordArray = 'None'
tableToPoints1Display.SelectNormalArray = 'None'
tableToPoints1Display.SelectTangentArray = 'None'
tableToPoints1Display.OSPRayScaleArray = 'TSL_norm'
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 0.19471644
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.GaussianRadius = 0.009735822
tableToPoints1Display.SetScaleArray = ['POINTS', 'TSL_norm']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = ['POINTS', 'TSL_norm']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'
tableToPoints1Display.SelectInputVectors = [None, '']
tableToPoints1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.00450234, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.00450234, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4999998, 0.9999977999999999, 10000.0]
renderView1.CameraFocalPoint = [0.4999998, 0.9999977999999999, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 8.0

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 6.0

# create a new 'XML Partitioned Unstructured Grid Reader'
solution_001pvtu = XMLPartitionedUnstructuredGridReader(registrationName='solution_001.pvtu*', FileName=['/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/output/solution_001.pvtu', '/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/output/solution_002.pvtu'])
solution_001pvtu.PointArrayStatus = ['u', 'strain', 'stress', 'subdomain']

# Properties modified on solution_001pvtu
solution_001pvtu.TimeArray = 'None'

# show data in view
solution_001pvtuDisplay = Show(solution_001pvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
solution_001pvtuDisplay.Representation = 'Surface'
solution_001pvtuDisplay.ColorArrayName = [None, '']
solution_001pvtuDisplay.SelectTCoordArray = 'None'
solution_001pvtuDisplay.SelectNormalArray = 'None'
solution_001pvtuDisplay.SelectTangentArray = 'None'
solution_001pvtuDisplay.OSPRayScaleArray = 'strain'
solution_001pvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
solution_001pvtuDisplay.SelectOrientationVectors = 'None'
solution_001pvtuDisplay.ScaleFactor = 0.2
solution_001pvtuDisplay.SelectScaleArray = 'None'
solution_001pvtuDisplay.GlyphType = 'Arrow'
solution_001pvtuDisplay.GlyphTableIndexArray = 'None'
solution_001pvtuDisplay.GaussianRadius = 0.01
solution_001pvtuDisplay.SetScaleArray = ['POINTS', 'strain']
solution_001pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution_001pvtuDisplay.OpacityArray = ['POINTS', 'strain']
solution_001pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
solution_001pvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
solution_001pvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
solution_001pvtuDisplay.ScalarOpacityUnitDistance = 0.44369208246944625
solution_001pvtuDisplay.OpacityArrayName = ['POINTS', 'strain']
solution_001pvtuDisplay.SelectInputVectors = ['POINTS', 'u']
solution_001pvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
solution_001pvtuDisplay.ScaleTransferFunction.Points = [-0.0004607821174431592, 0.0, 0.5, 0.0, 1.2898989325549337e-06, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
solution_001pvtuDisplay.OpacityTransferFunction.Points = [-0.0004607821174431592, 0.0, 0.5, 0.0, 1.2898989325549337e-06, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change representation type
solution_001pvtuDisplay.SetRepresentationType('Surface With Edges')

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(registrationName='WarpByVector1', Input=solution_001pvtu)
warpByVector1.Vectors = ['POINTS', 'u']

# show data in view
warpByVector1Display = Show(warpByVector1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.ColorArrayName = [None, '']
warpByVector1Display.SelectTCoordArray = 'None'
warpByVector1Display.SelectNormalArray = 'None'
warpByVector1Display.SelectTangentArray = 'None'
warpByVector1Display.OSPRayScaleArray = 'strain'
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.SelectOrientationVectors = 'None'
warpByVector1Display.ScaleFactor = 0.20019973290010284
warpByVector1Display.SelectScaleArray = 'None'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'None'
warpByVector1Display.GaussianRadius = 0.010009986645005143
warpByVector1Display.SetScaleArray = ['POINTS', 'strain']
warpByVector1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByVector1Display.OpacityArray = ['POINTS', 'strain']
warpByVector1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByVector1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByVector1Display.PolarAxes = 'PolarAxesRepresentation'
warpByVector1Display.ScalarOpacityUnitDistance = 0.4440463006157371
warpByVector1Display.OpacityArrayName = ['POINTS', 'strain']
warpByVector1Display.SelectInputVectors = ['POINTS', 'u']
warpByVector1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
warpByVector1Display.ScaleTransferFunction.Points = [-0.0004607821174431592, 0.0, 0.5, 0.0, 1.2898989325549337e-06, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
warpByVector1Display.OpacityTransferFunction.Points = [-0.0004607821174431592, 0.0, 0.5, 0.0, 1.2898989325549337e-06, 1.0, 0.5, 0.0]

# hide data in view
Hide(solution_001pvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(warpByVector1Display, ('POINTS', 'u', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
warpByVector1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1477, 789)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4999998, 0.9999977999999999, 10000.0]
renderView1.CameraFocalPoint = [0.4999998, 0.9999977999999999, 0.0]
renderView1.CameraParallelScale = 1.0826561294549992

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).