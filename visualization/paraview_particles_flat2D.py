try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

#############################################################

# PARAMETER visu
f =open('../bin/PARAMETER_MicroVic_flat.txt')
lines=f.readlines()
f.close()
Lx = float(lines[12-1])
Ly = float(lines[13-1])
dt = float(lines[18-1])
jump = float(lines[37-1])
dt = dt*jump

scaleU = .4

# view
RenderView1 = GetRenderView()
RenderView1.CenterAxesVisibility = 0
RenderView1.OrientationAxesVisibility = 0
RenderView1.Background = [.333333333333, .3333333333333, .5]

# axis: [Lx,Ly]
particle_0000 = FindSource("particle_000*")
DataRepresentation1 = GetDisplayProperties(particle_0000)
DataRepresentation1.CustomBoundsActive = [1, 1, 1]
DataRepresentation1.CustomBounds = [0.0, Lx, 0.0, Ly, 0.0, 0.0]
DataRepresentation1.CubeAxesVisibility = 1

# glyphe: scaleU
particle_00000 = GetActiveSource()
Glyph1 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )
Glyph1.Vectors = ['POINTS', 'vectors']
Glyph1.GlyphType = "2D Glyph"
Glyph1.GlyphTransform.Scale = [scaleU, scaleU, scaleU]
#Glyph1.RandomMode = 0

DataRepresentation1 = Show()
Render()


# time scale: dt
rho2D_0000 = GetActiveSource()
AnnotateTimeFilter1 = AnnotateTimeFilter()
AnnotateTimeFilter1.Scale = dt
AnnotateTimeFilter1.Format = 'Time: %5.2f'
DataRepresentation3 = Show()
DataRepresentation3.Position = [0.4, 0.93]
