# trace generated using paraview version 5.8.0-RC1
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

### My additions ###
import os
import sys
filepath = sys.argv[1]

if filepath[-1]!='/':
    filepath += '/'
### End of my additions ###

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
simul_Particlespvd = PVDReader(FileName=filepath+'simul_Particles.pvd')
simul_Particlespvd.CellArrays = ['NormU', 'NormOm', 'CoordNumb']
simul_Particlespvd.PointArrays = []
simul_Particlespvd.ColumnArrays = []

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1545, 803]

# get layout
layout1 = GetLayout()

# show data in view
simul_ParticlespvdDisplay = Show(simul_Particlespvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
simul_ParticlespvdDisplay.Representation = 'Surface'
simul_ParticlespvdDisplay.ColorArrayName = [None, '']
simul_ParticlespvdDisplay.LookupTable = None
simul_ParticlespvdDisplay.MapScalars = 1
simul_ParticlespvdDisplay.MultiComponentsMapping = 0
simul_ParticlespvdDisplay.InterpolateScalarsBeforeMapping = 1
simul_ParticlespvdDisplay.Opacity = 1.0
simul_ParticlespvdDisplay.PointSize = 2.0
simul_ParticlespvdDisplay.LineWidth = 1.0
simul_ParticlespvdDisplay.RenderLinesAsTubes = 0
simul_ParticlespvdDisplay.RenderPointsAsSpheres = 0
simul_ParticlespvdDisplay.Interpolation = 'Gouraud'
simul_ParticlespvdDisplay.Specular = 0.0
simul_ParticlespvdDisplay.SpecularColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.SpecularPower = 100.0
simul_ParticlespvdDisplay.Luminosity = 0.0
simul_ParticlespvdDisplay.Ambient = 0.0
simul_ParticlespvdDisplay.Diffuse = 1.0
simul_ParticlespvdDisplay.Roughness = 0.3
simul_ParticlespvdDisplay.Metallic = 0.0
simul_ParticlespvdDisplay.Texture = None
simul_ParticlespvdDisplay.RepeatTextures = 1
simul_ParticlespvdDisplay.InterpolateTextures = 0
simul_ParticlespvdDisplay.SeamlessU = 0
simul_ParticlespvdDisplay.SeamlessV = 0
simul_ParticlespvdDisplay.UseMipmapTextures = 0
simul_ParticlespvdDisplay.BaseColorTexture = None
simul_ParticlespvdDisplay.NormalTexture = None
simul_ParticlespvdDisplay.NormalScale = 1.0
simul_ParticlespvdDisplay.MaterialTexture = None
simul_ParticlespvdDisplay.OcclusionStrength = 1.0
simul_ParticlespvdDisplay.EmissiveTexture = None
simul_ParticlespvdDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.FlipTextures = 0
simul_ParticlespvdDisplay.BackfaceRepresentation = 'Follow Frontface'
simul_ParticlespvdDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.BackfaceOpacity = 1.0
simul_ParticlespvdDisplay.Position = [0.0, 0.0, 0.0]
simul_ParticlespvdDisplay.Scale = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.Orientation = [0.0, 0.0, 0.0]
simul_ParticlespvdDisplay.Origin = [0.0, 0.0, 0.0]
simul_ParticlespvdDisplay.Pickable = 1
simul_ParticlespvdDisplay.Triangulate = 0
simul_ParticlespvdDisplay.UseShaderReplacements = 0
simul_ParticlespvdDisplay.ShaderReplacements = ''
simul_ParticlespvdDisplay.NonlinearSubdivisionLevel = 1
simul_ParticlespvdDisplay.UseDataPartitions = 0
simul_ParticlespvdDisplay.OSPRayUseScaleArray = 0
simul_ParticlespvdDisplay.OSPRayScaleArray = ''
simul_ParticlespvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
simul_ParticlespvdDisplay.OSPRayMaterial = 'None'
simul_ParticlespvdDisplay.Orient = 0
simul_ParticlespvdDisplay.OrientationMode = 'Direction'
simul_ParticlespvdDisplay.SelectOrientationVectors = 'None'
simul_ParticlespvdDisplay.Scaling = 0
simul_ParticlespvdDisplay.ScaleMode = 'No Data Scaling Off'
simul_ParticlespvdDisplay.ScaleFactor = 0.02843133756687166
simul_ParticlespvdDisplay.SelectScaleArray = 'None'
simul_ParticlespvdDisplay.GlyphType = 'Arrow'
simul_ParticlespvdDisplay.UseGlyphTable = 0
simul_ParticlespvdDisplay.GlyphTableIndexArray = 'None'
simul_ParticlespvdDisplay.UseCompositeGlyphTable = 0
simul_ParticlespvdDisplay.UseGlyphCullingAndLOD = 0
simul_ParticlespvdDisplay.LODValues = []
simul_ParticlespvdDisplay.ColorByLODIndex = 0
simul_ParticlespvdDisplay.GaussianRadius = 0.0014215668783435832
simul_ParticlespvdDisplay.ShaderPreset = 'Sphere'
simul_ParticlespvdDisplay.CustomTriangleScale = 3
simul_ParticlespvdDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
simul_ParticlespvdDisplay.Emissive = 0
simul_ParticlespvdDisplay.ScaleByArray = 0
simul_ParticlespvdDisplay.SetScaleArray = [None, '']
simul_ParticlespvdDisplay.ScaleArrayComponent = 0
simul_ParticlespvdDisplay.UseScaleFunction = 1
simul_ParticlespvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
simul_ParticlespvdDisplay.OpacityByArray = 0
simul_ParticlespvdDisplay.OpacityArray = [None, '']
simul_ParticlespvdDisplay.OpacityArrayComponent = 0
simul_ParticlespvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
simul_ParticlespvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
simul_ParticlespvdDisplay.SelectionCellLabelBold = 0
simul_ParticlespvdDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
simul_ParticlespvdDisplay.SelectionCellLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.SelectionCellLabelFontFile = ''
simul_ParticlespvdDisplay.SelectionCellLabelFontSize = 18
simul_ParticlespvdDisplay.SelectionCellLabelItalic = 0
simul_ParticlespvdDisplay.SelectionCellLabelJustification = 'Left'
simul_ParticlespvdDisplay.SelectionCellLabelOpacity = 1.0
simul_ParticlespvdDisplay.SelectionCellLabelShadow = 0
simul_ParticlespvdDisplay.SelectionPointLabelBold = 0
simul_ParticlespvdDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
simul_ParticlespvdDisplay.SelectionPointLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.SelectionPointLabelFontFile = ''
simul_ParticlespvdDisplay.SelectionPointLabelFontSize = 18
simul_ParticlespvdDisplay.SelectionPointLabelItalic = 0
simul_ParticlespvdDisplay.SelectionPointLabelJustification = 'Left'
simul_ParticlespvdDisplay.SelectionPointLabelOpacity = 1.0
simul_ParticlespvdDisplay.SelectionPointLabelShadow = 0
simul_ParticlespvdDisplay.PolarAxes = 'PolarAxesRepresentation'
simul_ParticlespvdDisplay.ScalarOpacityFunction = None
simul_ParticlespvdDisplay.ScalarOpacityUnitDistance = 0.023751538017530917
simul_ParticlespvdDisplay.ExtractedBlockIndex = 0
simul_ParticlespvdDisplay.SelectMapper = 'Projected tetra'
simul_ParticlespvdDisplay.SamplingDimensions = [128, 128, 128]
simul_ParticlespvdDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
simul_ParticlespvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
simul_ParticlespvdDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
simul_ParticlespvdDisplay.GlyphType.TipResolution = 6
simul_ParticlespvdDisplay.GlyphType.TipRadius = 0.1
simul_ParticlespvdDisplay.GlyphType.TipLength = 0.35
simul_ParticlespvdDisplay.GlyphType.ShaftResolution = 6
simul_ParticlespvdDisplay.GlyphType.ShaftRadius = 0.03
simul_ParticlespvdDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
simul_ParticlespvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
simul_ParticlespvdDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
simul_ParticlespvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
simul_ParticlespvdDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
simul_ParticlespvdDisplay.DataAxesGrid.XTitle = 'X Axis'
simul_ParticlespvdDisplay.DataAxesGrid.YTitle = 'Y Axis'
simul_ParticlespvdDisplay.DataAxesGrid.ZTitle = 'Z Axis'
simul_ParticlespvdDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.XTitleFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.XTitleBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.XTitleItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.XTitleFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.XTitleShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.XTitleOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.YTitleFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.YTitleBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.YTitleItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.YTitleFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.YTitleShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.YTitleOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZTitleOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.FacesToRender = 63
simul_ParticlespvdDisplay.DataAxesGrid.CullBackface = 0
simul_ParticlespvdDisplay.DataAxesGrid.CullFrontface = 1
simul_ParticlespvdDisplay.DataAxesGrid.ShowGrid = 0
simul_ParticlespvdDisplay.DataAxesGrid.ShowEdges = 1
simul_ParticlespvdDisplay.DataAxesGrid.ShowTicks = 1
simul_ParticlespvdDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
simul_ParticlespvdDisplay.DataAxesGrid.AxesToLabel = 63
simul_ParticlespvdDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.XLabelFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.XLabelBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.XLabelItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.XLabelFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.XLabelShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.XLabelOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.YLabelFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.YLabelBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.YLabelItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.YLabelFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.YLabelShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.YLabelOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelFontFile = ''
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelBold = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelItalic = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelFontSize = 12
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelShadow = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZLabelOpacity = 1.0
simul_ParticlespvdDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
simul_ParticlespvdDisplay.DataAxesGrid.XAxisPrecision = 2
simul_ParticlespvdDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
simul_ParticlespvdDisplay.DataAxesGrid.XAxisLabels = []
simul_ParticlespvdDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
simul_ParticlespvdDisplay.DataAxesGrid.YAxisPrecision = 2
simul_ParticlespvdDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
simul_ParticlespvdDisplay.DataAxesGrid.YAxisLabels = []
simul_ParticlespvdDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
simul_ParticlespvdDisplay.DataAxesGrid.ZAxisPrecision = 2
simul_ParticlespvdDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
simul_ParticlespvdDisplay.DataAxesGrid.ZAxisLabels = []
simul_ParticlespvdDisplay.DataAxesGrid.UseCustomBounds = 0
simul_ParticlespvdDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
simul_ParticlespvdDisplay.PolarAxes.Visibility = 0
simul_ParticlespvdDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
simul_ParticlespvdDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
simul_ParticlespvdDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
simul_ParticlespvdDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.EnableCustomRange = 0
simul_ParticlespvdDisplay.PolarAxes.CustomRange = [0.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.PolarAxisVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.RadialAxesVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.DrawRadialGridlines = 1
simul_ParticlespvdDisplay.PolarAxes.PolarArcsVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.DrawPolarArcsGridlines = 1
simul_ParticlespvdDisplay.PolarAxes.NumberOfRadialAxes = 0
simul_ParticlespvdDisplay.PolarAxes.AutoSubdividePolarAxis = 1
simul_ParticlespvdDisplay.PolarAxes.NumberOfPolarAxis = 0
simul_ParticlespvdDisplay.PolarAxes.MinimumRadius = 0.0
simul_ParticlespvdDisplay.PolarAxes.MinimumAngle = 0.0
simul_ParticlespvdDisplay.PolarAxes.MaximumAngle = 90.0
simul_ParticlespvdDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
simul_ParticlespvdDisplay.PolarAxes.Ratio = 1.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
simul_ParticlespvdDisplay.PolarAxes.PolarLabelVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
simul_ParticlespvdDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
simul_ParticlespvdDisplay.PolarAxes.RadialLabelVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
simul_ParticlespvdDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
simul_ParticlespvdDisplay.PolarAxes.RadialUnitsVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.ScreenSize = 10.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleBold = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleItalic = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleShadow = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTitleFontSize = 12
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelBold = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelItalic = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelShadow = 0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisLabelFontSize = 12
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextBold = 0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextItalic = 0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextShadow = 0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
simul_ParticlespvdDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
simul_ParticlespvdDisplay.PolarAxes.EnableDistanceLOD = 1
simul_ParticlespvdDisplay.PolarAxes.DistanceLODThreshold = 0.7
simul_ParticlespvdDisplay.PolarAxes.EnableViewAngleLOD = 1
simul_ParticlespvdDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
simul_ParticlespvdDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
simul_ParticlespvdDisplay.PolarAxes.PolarTicksVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
simul_ParticlespvdDisplay.PolarAxes.TickLocation = 'Both'
simul_ParticlespvdDisplay.PolarAxes.AxisTickVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.AxisMinorTickVisibility = 0
simul_ParticlespvdDisplay.PolarAxes.ArcTickVisibility = 1
simul_ParticlespvdDisplay.PolarAxes.ArcMinorTickVisibility = 0
simul_ParticlespvdDisplay.PolarAxes.DeltaAngleMajor = 10.0
simul_ParticlespvdDisplay.PolarAxes.DeltaAngleMinor = 5.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
simul_ParticlespvdDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
simul_ParticlespvdDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
simul_ParticlespvdDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
simul_ParticlespvdDisplay.PolarAxes.ArcMajorTickSize = 0.0
simul_ParticlespvdDisplay.PolarAxes.ArcTickRatioSize = 0.3
simul_ParticlespvdDisplay.PolarAxes.ArcMajorTickThickness = 1.0
simul_ParticlespvdDisplay.PolarAxes.ArcTickRatioThickness = 0.5
simul_ParticlespvdDisplay.PolarAxes.Use2DMode = 0
simul_ParticlespvdDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

#change interaction mode for render view
renderView1.InteractionMode = '2D'

animationScene1.GoToLast()

# reset view to fit data bounds
renderView1.ResetCamera(-7.190294854808599e-05, 0.49449053406715393, -6.877483247080818e-05, 0.20007194578647614, -7.289636414498091e-05, 0.09390579164028168)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1
renderView1.Background = [0.0, 0.0, 0.0]

# Properties modified on simul_ParticlespvdDisplay
simul_ParticlespvdDisplay.Luminosity = 100.0

# Properties modified on simul_ParticlespvdDisplay
simul_ParticlespvdDisplay.Luminosity = 0.0

# Properties modified on simul_ParticlespvdDisplay
simul_ParticlespvdDisplay.Ambient = 1.

animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.24720931555930292, 0.10000158547700266, 1.093474335596475]
renderView1.CameraFocalPoint = [0.24720931555930292, 0.10000158547700266, 0.04691644763806835]
renderView1.CameraParallelScale = 0.1850072489624381

# save screenshot
SaveScreenshot(filepath+'initial_frame_top.png', renderView1, ImageResolution=[1545, 803],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    # PNG options
    CompressionLevel='5')

animationScene1.GoToLast()

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.24720931555930292, 0.10000158547700266, 1.093474335596475]
renderView1.CameraFocalPoint = [0.24720931555930292, 0.10000158547700266, 0.04691644763806835]
renderView1.CameraParallelScale = 0.1850072489624381

# save screenshot
SaveScreenshot(filepath+'last_frame_top.png', renderView1, ImageResolution=[1545, 803],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    # PNG options
    CompressionLevel='5')

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToFirst()

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.36931967543544053, -0.9465563024814039, 0.14339037836893895]
renderView1.CameraFocalPoint = [0.36931967543544053, 0.10000158547700266, 0.14339037836893895]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.2708691132059057

# save screenshot
SaveScreenshot(filepath+'initial_frame_side.png', renderView1, ImageResolution=[1545, 803],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    # PNG options
    CompressionLevel='5')

animationScene1.GoToLast()

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.36931967543544053, -0.9465563024814039, 0.14339037836893895]
renderView1.CameraFocalPoint = [0.36931967543544053, 0.10000158547700266, 0.14339037836893895]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.2708691132059057

# save screenshot
SaveScreenshot(filepath+'last_frame_side.png', renderView1, ImageResolution=[1545, 803],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    # PNG options
    CompressionLevel='5')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.36931967543544053, -0.9465563024814039, 0.14339037836893895]
renderView1.CameraFocalPoint = [0.36931967543544053, 0.10000158547700266, 0.14339037836893895]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.2708691132059057


#### SAVE WIDTH VIEW FROM HERE
# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.CameraPosition = [-0.5884626256622852, 0.0999999653067789, 0.12173497098410735]
renderView1.CameraFocalPoint = [0.0400000148001709, 0.0999999653067789, 0.12173497098410735]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.16265810048710172

# save screenshot
SaveScreenshot(filepath+'initial_frame_width.png', magnification=1, quality=100, view=renderView1)

animationScene1.GoToLast()

# current camera placement for renderView1
renderView1.CameraPosition = [-0.5884626256622852, 0.0999999653067789, 0.12173497098410735]
renderView1.CameraFocalPoint = [0.0400000148001709, 0.0999999653067789, 0.12173497098410735]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.16265810048710172

# save screenshot
SaveScreenshot(filepath+'last_frame_width.png', magnification=1, quality=100, view=renderView1)


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
