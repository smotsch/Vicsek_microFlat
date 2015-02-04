try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# PARAMETER visu
Lx = 10
Ly = 10
scaleU = 1
dt = .05


rho2D_00000 = GetActiveSource()
Tetrahedralize1 = Tetrahedralize()

RenderView1 = GetRenderView()
RenderView1.Background = [.333333333333, .3333333333333, .5]
RenderView1.InteractionMode = '3D'
DataRepresentation1 = GetDisplayProperties(rho2D_00000)
DataRepresentation1.Visibility = 0
a1_Density_PVLookupTable = GetLookupTableForArray( "Density", 1 )

DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = 'Density'
DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ColorArrayName = ('POINT_DATA', 'Density')
DataRepresentation2.ScalarOpacityUnitDistance = 0.9596915518332424
DataRepresentation2.LookupTable = a1_Density_PVLookupTable
DataRepresentation2.Visibility = 0

ExtractSurface1 = ExtractSurface()

DataRepresentation3 = Show()
DataRepresentation3.ColorArrayName = ('POINT_DATA', 'Density')
DataRepresentation3.LookupTable = a1_Density_PVLookupTable
DataRepresentation3.SelectionPointFieldDataArrayName = 'Density'
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

LoopSubdivision1 = LoopSubdivision()

DataRepresentation4 = Show()
DataRepresentation4.ColorArrayName = ('POINT_DATA', 'Density')
DataRepresentation4.LookupTable = a1_Density_PVLookupTable
DataRepresentation4.SelectionPointFieldDataArrayName = 'Density'
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]

DataRepresentation3.Visibility = 0

WarpByScalar1 = WarpByScalar()

WarpByScalar1.Scalars = ['POINTS', 'Density']

DataRepresentation5 = Show()
DataRepresentation5.ColorArrayName = ('POINT_DATA', 'Density')
DataRepresentation5.LookupTable = a1_Density_PVLookupTable
DataRepresentation5.SelectionPointFieldDataArrayName = 'Density'
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]

DataRepresentation4.Visibility = 0

WarpByScalar1.ScaleFactor = 1.0

GenerateSurfaceNormals1 = GenerateSurfaceNormals()

DataRepresentation6 = Show()
DataRepresentation6.ColorArrayName = ('POINT_DATA', 'Density')
DataRepresentation6.LookupTable = a1_Density_PVLookupTable
DataRepresentation6.SelectionPointFieldDataArrayName = 'Density'
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5000076295109483]

DataRepresentation5.Visibility = 0

# camera
# RenderView6 = GetRenderView()
# RenderView6.CameraPosition = [5.0, -12.0, 21.0]
# RenderView6.CameraViewAngle = 25.0
# RenderView6.CameraFocalPoint = [5.0, 5.0, 0.0]

# time scale: dt
rho2D_0000 = GetActiveSource()
AnnotateTimeFilter1 = AnnotateTimeFilter()
AnnotateTimeFilter1.Scale = dt
AnnotateTimeFilter1.Format = 'Time: %5.2f'
DataRepresentation7 = Show()
DataRepresentation7.Position = [0.4, 0.93]


##GenerateSurfaceNormals2 = FindSource("GenerateSurfaceNormals1")
my_representation4 = GetDisplayProperties(GenerateSurfaceNormals1)
my_representation4.CustomBoundsActive = [0, 0, 1]
my_representation4.CubeAxesZAxisTickVisibility = 0
my_representation4.CubeAxesZAxisVisibility = 0
my_representation4.CubeAxesZAxisMinorTickVisibility = 0
my_representation4.CubeAxesVisibility = 1
my_representation4.CustomRange = [0.0, 10.0, 0.0, 10.0, 0.102176, 0.732323]
my_representation4.CustomBounds = [0.0, 10.0, 0.0, 10.0, 0.0, 0.0]

Render()
