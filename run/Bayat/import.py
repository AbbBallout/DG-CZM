# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
solution0_0vtu = XMLUnstructuredGridReader(registrationName='solution0_0.vtu*', FileName=['/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/Bayat/solution0_0.vtu', '/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/Bayat/solution0_1.vtu'])
solution0_0vtu.PointArrayStatus = ['u', 'strain', 'stress']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on solution0_0vtu
solution0_0vtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
solution0_0vtuDisplay = Show(solution0_0vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
solution0_0vtuDisplay.Representation = 'Surface'
solution0_0vtuDisplay.ColorArrayName = [None, '']
solution0_0vtuDisplay.SelectTCoordArray = 'None'
solution0_0vtuDisplay.SelectNormalArray = 'None'
solution0_0vtuDisplay.SelectTangentArray = 'None'
solution0_0vtuDisplay.OSPRayScaleArray = 'strain'
solution0_0vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
solution0_0vtuDisplay.SelectOrientationVectors = 'None'
solution0_0vtuDisplay.ScaleFactor = 0.2
solution0_0vtuDisplay.SelectScaleArray = 'None'
solution0_0vtuDisplay.GlyphType = 'Arrow'
solution0_0vtuDisplay.GlyphTableIndexArray = 'None'
solution0_0vtuDisplay.GaussianRadius = 0.01
solution0_0vtuDisplay.SetScaleArray = ['POINTS', 'strain']
solution0_0vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution0_0vtuDisplay.OpacityArray = ['POINTS', 'strain']
solution0_0vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
solution0_0vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
solution0_0vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
solution0_0vtuDisplay.ScalarOpacityUnitDistance = 0.7043172784619601
solution0_0vtuDisplay.OpacityArrayName = ['POINTS', 'strain']
solution0_0vtuDisplay.SelectInputVectors = ['POINTS', 'u']
solution0_0vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
solution0_0vtuDisplay.ScaleTransferFunction.Points = [-0.00023086811415851116, 0.0, 0.5, 0.0, 1.2138886518187064e-07, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
solution0_0vtuDisplay.OpacityTransferFunction.Points = [-0.00023086811415851116, 0.0, 0.5, 0.0, 1.2138886518187064e-07, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.5, 1.0, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(solution0_0vtuDisplay, ('POINTS', 'u', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
solution0_0vtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
solution0_0vtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=solution0_0vtu)
plotOverLine1.Point2 = [1.0, 2.0, 0.0]

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [0.0, 1.0, 0.0]
plotOverLine1.Point2 = [1.0, 1.0, 0.0]

# show data in view
plotOverLine1Display = Show(plotOverLine1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plotOverLine1Display.Representation = 'Surface'
plotOverLine1Display.ColorArrayName = ['POINTS', 'u']
plotOverLine1Display.LookupTable = uLUT
plotOverLine1Display.SelectTCoordArray = 'None'
plotOverLine1Display.SelectNormalArray = 'None'
plotOverLine1Display.SelectTangentArray = 'None'
plotOverLine1Display.OSPRayScaleArray = 'arc_length'
plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine1Display.SelectOrientationVectors = 'None'
plotOverLine1Display.ScaleFactor = 0.1
plotOverLine1Display.SelectScaleArray = 'None'
plotOverLine1Display.GlyphType = 'Arrow'
plotOverLine1Display.GlyphTableIndexArray = 'None'
plotOverLine1Display.GaussianRadius = 0.005
plotOverLine1Display.SetScaleArray = ['POINTS', 'arc_length']
plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.OpacityArray = ['POINTS', 'arc_length']
plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'
plotOverLine1Display.SelectInputVectors = ['POINTS', 'u']
plotOverLine1Display.WriteLog = ''

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')

# trace defaults for the display properties.
plotOverLine1Display_1.UseIndexForXAxis = 0
plotOverLine1Display_1.XArrayName = 'arc_length'
plotOverLine1Display_1.SeriesVisibility = ['strain_Magnitude', 'stress_Magnitude', 'u_Magnitude']
plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'strain_0', 'strain_0', 'strain_1', 'strain_1', 'strain_2', 'strain_2', 'strain_3', 'strain_3', 'strain_4', 'strain_4', 'strain_5', 'strain_5', 'strain_6', 'strain_6', 'strain_7', 'strain_7', 'strain_8', 'strain_8', 'strain_Magnitude', 'strain_Magnitude', 'stress_0', 'stress_0', 'stress_1', 'stress_1', 'stress_2', 'stress_2', 'stress_3', 'stress_3', 'stress_4', 'stress_4', 'stress_5', 'stress_5', 'stress_6', 'stress_6', 'stress_7', 'stress_7', 'stress_8', 'stress_8', 'stress_Magnitude', 'stress_Magnitude', 'u_X', 'u_X', 'u_Y', 'u_Y', 'u_Z', 'u_Z', 'u_Magnitude', 'u_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'strain_0', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_1', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_2', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'strain_3', '0.6', '0.3100022888532845', '0.6399938963912413', 'strain_4', '1', '0.5000076295109483', '0', 'strain_5', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'strain_6', '0', '0', '0', 'strain_7', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_8', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_0', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_1', '1', '0.5000076295109483', '0', 'stress_2', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'stress_3', '0', '0', '0', 'stress_4', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'stress_5', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'stress_6', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_7', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_8', '1', '0.5000076295109483', '0', 'stress_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'u_X', '0', '0', '0', 'u_Y', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'u_Z', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'u_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'vtkValidPointMask', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_X', '1', '0.5000076295109483', '0', 'Points_Y', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'Points_Z', '0', '0', '0', 'Points_Magnitude', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845']
plotOverLine1Display_1.SeriesOpacity = ['arc_length', '1.0', 'strain_0', '1.0', 'strain_1', '1.0', 'strain_2', '1.0', 'strain_3', '1.0', 'strain_4', '1.0', 'strain_5', '1.0', 'strain_6', '1.0', 'strain_7', '1.0', 'strain_8', '1.0', 'strain_Magnitude', '1.0', 'stress_0', '1.0', 'stress_1', '1.0', 'stress_2', '1.0', 'stress_3', '1.0', 'stress_4', '1.0', 'stress_5', '1.0', 'stress_6', '1.0', 'stress_7', '1.0', 'stress_8', '1.0', 'stress_Magnitude', '1.0', 'u_X', '1.0', 'u_Y', '1.0', 'u_Z', '1.0', 'u_Magnitude', '1.0', 'vtkValidPointMask', '1.0', 'Points_X', '1.0', 'Points_Y', '1.0', 'Points_Z', '1.0', 'Points_Magnitude', '1.0']
plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesLabelPrefix = ''
plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'u_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesMarkerSize = ['arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'u_Magnitude', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesOpacity = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'vtkValidPointMask', '1']
plotOverLine1Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'vtkValidPointMask', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'vtkValidPointMask', '4']

# update the view to ensure updated data information
lineChartView1.Update()

# set active source
SetActiveSource(solution0_0vtu)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=plotOverLine1)

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(registrationName='PlotOverLine2', Input=solution0_0vtu)
plotOverLine2.Point2 = [1.0, 2.0, 0.0]

# set active view
SetActiveView(renderView1)

# Properties modified on plotOverLine2
plotOverLine2.Point1 = [0.0, 2.0, 0.0]

# show data in view
plotOverLine2Display = Show(plotOverLine2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plotOverLine2Display.Representation = 'Surface'
plotOverLine2Display.ColorArrayName = ['POINTS', 'u']
plotOverLine2Display.LookupTable = uLUT
plotOverLine2Display.SelectTCoordArray = 'None'
plotOverLine2Display.SelectNormalArray = 'None'
plotOverLine2Display.SelectTangentArray = 'None'
plotOverLine2Display.OSPRayScaleArray = 'arc_length'
plotOverLine2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine2Display.SelectOrientationVectors = 'None'
plotOverLine2Display.ScaleFactor = 0.1
plotOverLine2Display.SelectScaleArray = 'None'
plotOverLine2Display.GlyphType = 'Arrow'
plotOverLine2Display.GlyphTableIndexArray = 'None'
plotOverLine2Display.GaussianRadius = 0.005
plotOverLine2Display.SetScaleArray = ['POINTS', 'arc_length']
plotOverLine2Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine2Display.OpacityArray = ['POINTS', 'arc_length']
plotOverLine2Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine2Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine2Display.PolarAxes = 'PolarAxesRepresentation'
plotOverLine2Display.SelectInputVectors = ['POINTS', 'u']
plotOverLine2Display.WriteLog = ''

# Create a new 'Line Chart View'
lineChartView2 = CreateView('XYChartView')

# show data in view
plotOverLine2Display_1 = Show(plotOverLine2, lineChartView2, 'XYChartRepresentation')

# trace defaults for the display properties.
plotOverLine2Display_1.UseIndexForXAxis = 0
plotOverLine2Display_1.XArrayName = 'arc_length'
plotOverLine2Display_1.SeriesVisibility = ['strain_Magnitude', 'stress_Magnitude', 'u_Magnitude']
plotOverLine2Display_1.SeriesLabel = ['arc_length', 'arc_length', 'strain_0', 'strain_0', 'strain_1', 'strain_1', 'strain_2', 'strain_2', 'strain_3', 'strain_3', 'strain_4', 'strain_4', 'strain_5', 'strain_5', 'strain_6', 'strain_6', 'strain_7', 'strain_7', 'strain_8', 'strain_8', 'strain_Magnitude', 'strain_Magnitude', 'stress_0', 'stress_0', 'stress_1', 'stress_1', 'stress_2', 'stress_2', 'stress_3', 'stress_3', 'stress_4', 'stress_4', 'stress_5', 'stress_5', 'stress_6', 'stress_6', 'stress_7', 'stress_7', 'stress_8', 'stress_8', 'stress_Magnitude', 'stress_Magnitude', 'u_X', 'u_X', 'u_Y', 'u_Y', 'u_Z', 'u_Z', 'u_Magnitude', 'u_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine2Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'strain_0', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_1', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_2', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'strain_3', '0.6', '0.3100022888532845', '0.6399938963912413', 'strain_4', '1', '0.5000076295109483', '0', 'strain_5', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'strain_6', '0', '0', '0', 'strain_7', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_8', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_0', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_1', '1', '0.5000076295109483', '0', 'stress_2', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'stress_3', '0', '0', '0', 'stress_4', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'stress_5', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'stress_6', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_7', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_8', '1', '0.5000076295109483', '0', 'stress_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'u_X', '0', '0', '0', 'u_Y', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'u_Z', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'u_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'vtkValidPointMask', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_X', '1', '0.5000076295109483', '0', 'Points_Y', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'Points_Z', '0', '0', '0', 'Points_Magnitude', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845']
plotOverLine2Display_1.SeriesOpacity = ['arc_length', '1.0', 'strain_0', '1.0', 'strain_1', '1.0', 'strain_2', '1.0', 'strain_3', '1.0', 'strain_4', '1.0', 'strain_5', '1.0', 'strain_6', '1.0', 'strain_7', '1.0', 'strain_8', '1.0', 'strain_Magnitude', '1.0', 'stress_0', '1.0', 'stress_1', '1.0', 'stress_2', '1.0', 'stress_3', '1.0', 'stress_4', '1.0', 'stress_5', '1.0', 'stress_6', '1.0', 'stress_7', '1.0', 'stress_8', '1.0', 'stress_Magnitude', '1.0', 'u_X', '1.0', 'u_Y', '1.0', 'u_Z', '1.0', 'u_Magnitude', '1.0', 'vtkValidPointMask', '1.0', 'Points_X', '1.0', 'Points_Y', '1.0', 'Points_Z', '1.0', 'Points_Magnitude', '1.0']
plotOverLine2Display_1.SeriesPlotCorner = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesLabelPrefix = ''
plotOverLine2Display_1.SeriesLineStyle = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesLineThickness = ['arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'u_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesMarkerSize = ['arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'u_Magnitude', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView2, layout=layout1, hint=1)

# Properties modified on plotOverLine2Display_1
plotOverLine2Display_1.SeriesOpacity = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine2Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'vtkValidPointMask', '1']
plotOverLine2Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'vtkValidPointMask', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine2Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'vtkValidPointMask', '4']

# update the view to ensure updated data information
lineChartView1.Update()

# update the view to ensure updated data information
lineChartView2.Update()

# destroy lineChartView2
Delete(lineChartView2)
del lineChartView2

# close an empty frame
layout1.Collapse(4)

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(lineChartView1)

# set active source
SetActiveSource(plotOverLine2)

# show data in view
plotOverLine2Display_1 = Show(plotOverLine2, lineChartView1, 'XYChartRepresentation')

# trace defaults for the display properties.
plotOverLine2Display_1.UseIndexForXAxis = 0
plotOverLine2Display_1.XArrayName = 'arc_length'
plotOverLine2Display_1.SeriesVisibility = ['strain_Magnitude', 'stress_Magnitude', 'u_Magnitude']
plotOverLine2Display_1.SeriesLabel = ['arc_length', 'arc_length', 'strain_0', 'strain_0', 'strain_1', 'strain_1', 'strain_2', 'strain_2', 'strain_3', 'strain_3', 'strain_4', 'strain_4', 'strain_5', 'strain_5', 'strain_6', 'strain_6', 'strain_7', 'strain_7', 'strain_8', 'strain_8', 'strain_Magnitude', 'strain_Magnitude', 'stress_0', 'stress_0', 'stress_1', 'stress_1', 'stress_2', 'stress_2', 'stress_3', 'stress_3', 'stress_4', 'stress_4', 'stress_5', 'stress_5', 'stress_6', 'stress_6', 'stress_7', 'stress_7', 'stress_8', 'stress_8', 'stress_Magnitude', 'stress_Magnitude', 'u_X', 'u_X', 'u_Y', 'u_Y', 'u_Z', 'u_Z', 'u_Magnitude', 'u_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine2Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'strain_0', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_1', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_2', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'strain_3', '0.6', '0.3100022888532845', '0.6399938963912413', 'strain_4', '1', '0.5000076295109483', '0', 'strain_5', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'strain_6', '0', '0', '0', 'strain_7', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'strain_8', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'strain_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_0', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_1', '1', '0.5000076295109483', '0', 'stress_2', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'stress_3', '0', '0', '0', 'stress_4', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'stress_5', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'stress_6', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'stress_7', '0.6', '0.3100022888532845', '0.6399938963912413', 'stress_8', '1', '0.5000076295109483', '0', 'stress_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'u_X', '0', '0', '0', 'u_Y', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'u_Z', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'u_Magnitude', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'vtkValidPointMask', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_X', '1', '0.5000076295109483', '0', 'Points_Y', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'Points_Z', '0', '0', '0', 'Points_Magnitude', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845']
plotOverLine2Display_1.SeriesOpacity = ['arc_length', '1.0', 'strain_0', '1.0', 'strain_1', '1.0', 'strain_2', '1.0', 'strain_3', '1.0', 'strain_4', '1.0', 'strain_5', '1.0', 'strain_6', '1.0', 'strain_7', '1.0', 'strain_8', '1.0', 'strain_Magnitude', '1.0', 'stress_0', '1.0', 'stress_1', '1.0', 'stress_2', '1.0', 'stress_3', '1.0', 'stress_4', '1.0', 'stress_5', '1.0', 'stress_6', '1.0', 'stress_7', '1.0', 'stress_8', '1.0', 'stress_Magnitude', '1.0', 'u_X', '1.0', 'u_Y', '1.0', 'u_Z', '1.0', 'u_Magnitude', '1.0', 'vtkValidPointMask', '1.0', 'Points_X', '1.0', 'Points_Y', '1.0', 'Points_Z', '1.0', 'Points_Magnitude', '1.0']
plotOverLine2Display_1.SeriesPlotCorner = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesLabelPrefix = ''
plotOverLine2Display_1.SeriesLineStyle = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesLineThickness = ['arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'u_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'u_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesMarkerSize = ['arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'u_Magnitude', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

# Properties modified on plotOverLine2Display_1
plotOverLine2Display_1.SeriesOpacity = ['arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'u_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine2Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'arc_length', '1', 'strain_0', '1', 'strain_1', '1', 'strain_2', '1', 'strain_3', '1', 'strain_4', '1', 'strain_5', '1', 'strain_6', '1', 'strain_7', '1', 'strain_8', '1', 'strain_Magnitude', '1', 'stress_0', '1', 'stress_1', '1', 'stress_2', '1', 'stress_3', '1', 'stress_4', '1', 'stress_5', '1', 'stress_6', '1', 'stress_7', '1', 'stress_8', '1', 'stress_Magnitude', '1', 'u_Magnitude', '1', 'u_X', '1', 'u_Y', '1', 'u_Z', '1', 'vtkValidPointMask', '1']
plotOverLine2Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'arc_length', '2', 'strain_0', '2', 'strain_1', '2', 'strain_2', '2', 'strain_3', '2', 'strain_4', '2', 'strain_5', '2', 'strain_6', '2', 'strain_7', '2', 'strain_8', '2', 'strain_Magnitude', '2', 'stress_0', '2', 'stress_1', '2', 'stress_2', '2', 'stress_3', '2', 'stress_4', '2', 'stress_5', '2', 'stress_6', '2', 'stress_7', '2', 'stress_8', '2', 'stress_Magnitude', '2', 'u_Magnitude', '2', 'u_X', '2', 'u_Y', '2', 'u_Z', '2', 'vtkValidPointMask', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'arc_length', '0', 'strain_0', '0', 'strain_1', '0', 'strain_2', '0', 'strain_3', '0', 'strain_4', '0', 'strain_5', '0', 'strain_6', '0', 'strain_7', '0', 'strain_8', '0', 'strain_Magnitude', '0', 'stress_0', '0', 'stress_1', '0', 'stress_2', '0', 'stress_3', '0', 'stress_4', '0', 'stress_5', '0', 'stress_6', '0', 'stress_7', '0', 'stress_8', '0', 'stress_Magnitude', '0', 'u_Magnitude', '0', 'u_X', '0', 'u_Y', '0', 'u_Z', '0', 'vtkValidPointMask', '0']
plotOverLine2Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'arc_length', '4', 'strain_0', '4', 'strain_1', '4', 'strain_2', '4', 'strain_3', '4', 'strain_4', '4', 'strain_5', '4', 'strain_6', '4', 'strain_7', '4', 'strain_8', '4', 'strain_Magnitude', '4', 'stress_0', '4', 'stress_1', '4', 'stress_2', '4', 'stress_3', '4', 'stress_4', '4', 'stress_5', '4', 'stress_6', '4', 'stress_7', '4', 'stress_8', '4', 'stress_Magnitude', '4', 'u_Magnitude', '4', 'u_X', '4', 'u_Y', '4', 'u_Z', '4', 'vtkValidPointMask', '4']

# set active source
SetActiveSource(plotOverLine2)

# set active source
SetActiveSource(solution0_0vtu)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=plotOverLine2)

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(registrationName='WarpByVector1', Input=solution0_0vtu)
warpByVector1.Vectors = ['POINTS', 'u']

# set active view
SetActiveView(renderView1)

# hide data in view
Hide(solution0_0vtu, renderView1)

# set active source
SetActiveSource(warpByVector1)

# show data in view
warpByVector1Display = Show(warpByVector1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.ColorArrayName = ['POINTS', 'u']
warpByVector1Display.LookupTable = uLUT
warpByVector1Display.SelectTCoordArray = 'None'
warpByVector1Display.SelectNormalArray = 'None'
warpByVector1Display.SelectTangentArray = 'None'
warpByVector1Display.OSPRayScaleArray = 'strain'
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.SelectOrientationVectors = 'None'
warpByVector1Display.ScaleFactor = 0.2000999910860342
warpByVector1Display.SelectScaleArray = 'None'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'None'
warpByVector1Display.GaussianRadius = 0.010004999554301711
warpByVector1Display.SetScaleArray = ['POINTS', 'strain']
warpByVector1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByVector1Display.OpacityArray = ['POINTS', 'strain']
warpByVector1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByVector1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByVector1Display.PolarAxes = 'PolarAxesRepresentation'
warpByVector1Display.ScalarOpacityFunction = uPWF
warpByVector1Display.ScalarOpacityUnitDistance = 0.7045989919311844
warpByVector1Display.OpacityArrayName = ['POINTS', 'strain']
warpByVector1Display.SelectInputVectors = ['POINTS', 'u']
warpByVector1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
warpByVector1Display.ScaleTransferFunction.Points = [-0.00023086811415851116, 0.0, 0.5, 0.0, 1.2138886518187064e-07, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
warpByVector1Display.OpacityTransferFunction.Points = [-0.00023086811415851116, 0.0, 0.5, 0.0, 1.2138886518187064e-07, 1.0, 0.5, 0.0]

# show color bar/color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1529, 789)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 1.0, 10000.0]
renderView1.CameraFocalPoint = [0.5, 1.0, 0.0]
renderView1.CameraParallelScale = 1.118033988749895

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).