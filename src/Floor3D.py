########################################################################################################
### TOPOLOGY OPTIMIZATION FOR SIMPLY 3D FLOOR STRUCTURE                                              ###
###                                                                                                  ###
### Arjen Deetman                                                                                    ###
### version 17-09-2018                                                                               ###
########################################################################################################

########################################################################################################
### INPUT PARAMETERS                                                                                 ###
########################################################################################################

# File and job names
jobName='MODEL-CHECK' ## Job name to check the model
taskName='OPTIMIZATION-TASK' ## Task name for topology optimization
processName='OPT-3D-MODEL' # Process name for topology optimization
optJobName='{}-JOB'.format(processName) ## Job name for topology optimization

# Mesh size
meshSize=float(20.0) ## mm

# Element type
elType='C3D8' ## C3D8, C3D8R, C3D20 or C3D20R

# Elements
nelx=int(160)
nely=int(160) # nely=1 for a 2D beam! 
nelz=int(5)

# Geometry
sizeBC=float(meshSize*10) ## Size of the surface of the BCs in mm
fixy=False ## Fix all nodes in y-direction (True or False)

# Geometry restrictions
rmin=float(1.5*meshSize) ## Filter radius
minSizeTop=float(0.0) ## Minimum number of elements of the top layer
freezeCellsBC=True ## Freeze cells at boundary condition
planarSym1=True ## Planar symmetry for axis 1
planarSym2=True ## Planar symmetry for axis 2
demold=True ## Demold control
demoldAngle=float(0.0) ## Draft angle for demolding in degrees (0-20)

# Material properties
E=float(30000) ## N/mm2
v=float(0.2) ## -

# Optimization parameters
volumeFraction=float(0.5) ## Volume fraction range 0..1 (for design space)
operation='SUM' ## MAX, MIN or SUM: Operation between steps and load cases
maxStress=25 ## Maximum Von Mises stress (None is value means no stress contraint) 

# Computer settings
numCores=int(2) # Number of cores that will be uese of the processor
perMem=int(90) # Percent of the memory that will be uses of the computer

# Optimization task settings
## Material interpolation penalty
interpolationPenalty=float(3.0) ## Material interpolation penalty (default=3.0)
## Density
minDens=float(0.001) ## Minimum density (default=0.001)
maxDens=float(1.0) ## Maximum density (default=0.001)
maxDeltaDens=float(0.2) ## Maximum change per design cycle (default=0.25)
initialDens=volumeFraction ## Initial density (default=DEFAULT)
## Convergence
### Stops when both criteria are reached
conObjective=float(0.001) ## Objective function delta criterion (default=0.001)
conDensity=float(0.005) ## Element density delta criterion (default=0.005)
## Cycles
mDC=int(999) ## Maximum number of cycles

# Submit
submitJob=False ## Normal job to check the model (only for LC-1)
submitOpt=False ## Submitting the optimization task


########################################################################################################
### DEFINE LOADCASES AND LOADS                                                                       ###
########################################################################################################

# Number of load cases
nLoadCases=int(1)

# Number of tiles in top layer
nTilesX=int(2) ## Number of tiles in the top layer in x-direction
nTilesY=int(2) ## Number of tiles in the top layer in y-direction

# Set load for all tiles in top layer
for lc in range(1, nLoadCases+1):
    for tile in range(1, (nTilesX*nTilesY)+1):
        globals()["LC{}_T{}".format(lc, tile)]=0.0

# Set load for specific tile and load case
## Load case 1
LC1_T1=float(3000) ## Total force on tile
LC1_T2=float(3000) ## Total force on tile
LC1_T3=float(3000) ## Total force on tile
LC1_T4=float(3000) ## Total force on tile


########################################################################################################
### CREATE PART, MATERIALS, SECTION AND ASSEMBLY                                                     ###
########################################################################################################

# Loading modules
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

# First basic calculations
nElements=int(nelx*nely*nelz)
spanX=float(nelx*meshSize)
spanY=float(nely*meshSize)
height=float(nelz*meshSize)
if minSizeTop>0: 
    freezeTop=True
else: 
    freezeTop=False
minSizeTop=max(minSizeTop*meshSize, meshSize)

# Change tiles in y-direction
if nely==1: nTilesY=1

# Path
import os
folder=os.getcwd()
pathProcess=os.path.join(folder, processName)

# Always start with
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# Model name
modelName='LC-1' ## Name of the model (don't use Model-1)

# Create new model and delete Model-1
mdb.Model(modelType=STANDARD_EXPLICIT, name=modelName)
if 'Model-1' in mdb.models: del mdb.models['Model-1']

# Cleaning script
mod=mdb.models[modelName]

# Create Part
mod.ConstrainedSketch(name='__profile__', sheetSize=200.0)
## Cleaning script
modske=mod.sketches['__profile__']
## Line 1
p1=(0.0, 0.0)
p2=(spanX, 0.0)
modske.Line(point1=p1, point2=p2)
## Line 2
p1=p2
p2=(spanX, spanY)
modske.Line(point1=p1, point2=p2)
## Line 3
p1=p2
p2=(0.0, spanY)
modske.Line(point1=p1, point2=p2)
## Line 4
p1=p2
p2=(0.0, 0.0)
modske.Line(point1=p1, point2=p2)
## Create
mod.Part(dimensionality=THREE_D, name='PART-1', type=DEFORMABLE_BODY)
mod.parts['PART-1'].BaseSolidExtrude(depth=height, sketch=modske)
del modske, p1, p2

# Cleaning script
modpa=mod.parts['PART-1']

# Create sets
modpa.Set(cells=modpa.cells.getSequenceFromMask(('[#1 ]', ), ), name='PART-1')

# Cleaning script
modma=mod.materials

# Create material
mod.Material(name='CONCRETE')
modma['CONCRETE'].Elastic(table=((E, v), ))

# Create section
mod.HomogeneousSolidSection(material='CONCRETE', name='SECTION-1')

# Assign section
modpa.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
    region=modpa.sets['PART-1'], sectionName='SECTION-1', thicknessAssignment=FROM_SECTION)

# Cleaning script
modra=mod.rootAssembly

# Assembly
modra.DatumCsysByDefault(CARTESIAN)
modra.Instance(dependent=ON, name='PART-1-1', part=modpa)

# Cleaning script
modrain=modra.instances['PART-1-1']

# Create sets
# DESIGN-SPACE
findat=(spanX/2, spanY/2, height/2)
modra.Set(cells=modrain.cells.findAt((findat, )), name='DESIGN-SPACE')
## EDGE-1-TOP
findat=(0.0, spanY/2, height)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-1-TOP')
## EDGE-1-BOT
findat=(0.0, spanY/2, 0.0)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-1-BOT')
## EDGE-2-TOP
findat=(spanX, spanY/2, height)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-2-TOP')
## EDGE-2-BOT
findat=(spanX, spanY/2, 0.0)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-2-BOT')
## EDGE-3-TOP
findat=(spanX/2, 0.0, height)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-3-TOP')
## EDGE-3-BOT
findat=(spanX/2, 0.0, 0.0)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-3-BOT')
## EDGE-4-TOP
findat=(spanX/2, spanY, height)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-4-TOP')
## EDGE-4-BOT
findat=(spanX/2, spanY, 0.0)
modra.Set(edges=modrain.edges.findAt((findat, )), name='EDGE-4-BOT')
## SURF-TOP
findat=(spanX/2, spanY/2, height)
modra.Surface(name='SURF-TOP', side1Faces=modrain.faces.findAt((findat, )))
## SURF-BOT
findat=(spanX/2, spanY/2, 0.0)
modra.Surface(name='SURF-BOT', side1Faces=modrain.faces.findAt((findat, )))
del findat


########################################################################################################
### CREATE SETS FOR BC IN Y-DIRECTION                                                                ###
########################################################################################################

# Create sets if nodes in y-direction are fixed
if fixy==True:
    ## FACE-Y-1
    findat=(spanX/2, 0.0, height/2)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-Y-1')
    ## FACE-Y-2
    findat=(spanX/2, spanY, height/2)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-Y-2')


########################################################################################################
### FUNCTIONS FOR PARTIION OF EDGES AND CELLS                                                        ###
########################################################################################################

# Define function for partition of edge and cell
## Partition edge
def partitionEdge(fPoint):
    modpa.PartitionEdgeByPoint(edge=modpa.edges.findAt(fPoint, ), point=fPoint)
## Partition cell
def partitionCell(fCell, fEdge, fPoint):
    modpa.PartitionCellByPlanePointNormal(cells=modpa.cells.findAt((fCell, )), 
        normal=modpa.edges.findAt(fEdge, ), point=modpa.vertices.findAt(fPoint, ))


########################################################################################################
### SPLIT PART IN TOP-LAYER AND BOTTOM-LAYER                                                         ###
########################################################################################################

# Creating partition: split in bottom- and top-layer
## Define cell, edge and point
fCell=(spanX/2, spanY/2, height/2)
fPoint=(0.0, 0.0, height-minSizeTop)
fEdge=(0.0, 0.0, height/2)
## Apply functions
partitionEdge(fPoint)
partitionCell(fCell, fEdge, fPoint)
del fCell, fPoint, fEdge

# Create sets on assembly after splitting in top and bottom layer
## TOP-LAYER
findat=(spanX/2, spanY/2, height-float(minSizeTop)/2)
modpa.Set(cells=modpa.cells.findAt((findat, )), name='TOP-LAYER')
## BOTTOM-LAYER
findat=(spanX/2, spanY/2, float(height-minSizeTop)/2)
modpa.Set(cells=modpa.cells.findAt((findat, )), name='BOTTOM-LAYER')
## NODE-DEMOLD (node in demold plane)
findat=(0.0, 0.0, height)
modra.Set(name='NODE-DEMOLD', vertices=modrain.vertices.findAt((findat, )))
del findat


########################################################################################################
### CREATING TILES IN TOP-LAYER (surfaces for loads)                                                 ###
########################################################################################################

# Create tiles in top layer
## X-DIRECTION
if minSizeTop>0:
    for x in range(1, nTilesX):
        ## Define cell, edge and point
        fCell=(float(x-0.5)/nTilesX*spanX, spanY/2, height-float(minSizeTop)/2)
        fPoint=(float(x)/nTilesX*spanX, 0.0, height)
        fEdge=(float(x+0.25)/nTilesX*spanX, 0.0, height)
        ## Apply functions
        partitionEdge(fPoint)
        partitionCell(fCell, fEdge, fPoint)
        del fCell, fPoint, fEdge
## Y-DIRECTION
if minSizeTop>0:
    for x in range(0, nTilesX):
        for y in range(1, nTilesY):
            ## Define cell, edge and point
            fCell=(float(x+0.5)/nTilesX*spanX, float(y-0.5)/nTilesY*spanY, height-float(minSizeTop)/2)
            fPoint=(float(x+1)/nTilesX*spanX, float(y)/nTilesY*spanY, height)
            fEdge=(float(x+1)/nTilesX*spanX, float(y-0.25)/nTilesY*spanY, height)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge

# Create load surfaces
for y in range(0, nTilesY):
    for x in range(0, nTilesX):
        ## Define center of tile
        findat=((0.5+x)*spanX/nTilesX, (0.5+y)*spanY/nTilesY, height)
        n=(x+1)+(nTilesX*y)
        ## Create surface
        modpa.Surface(name='LOAD-TILE-{}'.format(n), side1Faces=modpa.faces.findAt((findat, )))
        del findat, n


########################################################################################################
### PARTITION OF DESIGN SPACE FOR BOUNDARY CONDITIONS                                                ###
########################################################################################################

# Divide design space in nine cells
# Divide bottom layer
if sizeBC!=0:
    ## X-DIRECTION
    for x in [sizeBC, spanX/2, spanX-sizeBC]:
        ## Define cell, edge and point
        fCell=(x-0.1, spanY/2, (height-minSizeTop)/2)
        fPoint=(x, 0.0, 0.0)
        fEdge=(x, 0.0, 0.0)
        ## Apply functions
        partitionEdge(fPoint)
        partitionCell(fCell, fEdge, fPoint)
        del fCell, fPoint, fEdge
if nely!=1 and sizeBC!=0:
    ## Y-DIRECTION
    for x in [sizeBC, spanX/2, spanX-sizeBC, spanX]:
        for y in [sizeBC, spanY/2, spanY-sizeBC]:
            ## Define cell, edge and point
            fCell=(x-0.1, y-0.1, (height-minSizeTop)/2)
            fPoint=(x, y, 0.0)
            fEdge=(x, y, 0.0)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge
if sizeBC==0:
    ## X-DIRECTION
    for x in [spanX/2]:
        ## Define cell, edge and point
        fCell=(x-0.1, spanY/2, (height-minSizeTop)/2)
        fPoint=(x, 0.0, 0.0)
        fEdge=(x, 0.0, 0.0)
        ## Apply functions
        partitionEdge(fPoint)
        partitionCell(fCell, fEdge, fPoint)
        del fCell, fPoint, fEdge
if sizeBC==0 and nely!=1:
    ## Y-DIRECTION
    for x in [spanX/2, spanX]:
        for y in [spanY/2]:
            ## Define cell, edge and point
            fCell=(x-0.1, y-0.1, (height-minSizeTop)/2)
            fPoint=(x, y, 0.0)
            fEdge=(x, y, 0.0)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge

# Create sets on assembly 
## NODE-EDGE-3-MID-TOP
findat=(spanX/2, 0.0, height)
modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-3-MID-TOP')
## NODE-EDGE-3-MID-BOT
findat=(spanX/2, 0.0, 0.0)
modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-3-MID-BOT')
## NODE-EDGE-4-MID-TOP
findat=(spanX/2, spanY, height)
modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-4-MID-TOP')
## NODE-EDGE-4-MID-BOT
findat=(spanX/2, spanY, 0.0)
modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-4-MID-BOT')
if nely==1 and sizeBC!=0:
    ## FACE-BC-1
    findat=(sizeBC/2, meshSize/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-1')
    ## FACE-BC-2
    findat=(spanX-sizeBC/2, meshSize/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-2')
if nely!=1:
    ## NODE-MIDDLE-TOP
    findat=(spanX/2, spanY/2, height)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-MIDDLE-TOP')
    ## NODE-MIDDLE-BOTTOM
    findat=(spanX/2, spanY/2, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-MIDDLE-BOT')
    ## NODE-EDGE-1-MID-TOP
    findat=(0.0, spanY/2, height)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-1-MID-TOP')
    ## NODE-EDGE-1-MID-BOT
    findat=(0.0, spanY/2, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-1-MID-BOT')
    ## NODE-EDGE-2-MID-TOP
    findat=(spanX, spanY/2, height)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-2-MID-TOP')
    ## NODE-EDGE-2-MID-BOT
    findat=(spanX, spanY/2, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-EDGE-2-MID-BOT')
if nely!=1 and sizeBC!=0:
    ## FACE-BC-1
    findat=(sizeBC/2, sizeBC/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-1')
    ## FACE-BC-2
    findat=(spanX-sizeBC/2, sizeBC/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-2')
    ## FACE-BC-3
    findat=(sizeBC/2, spanY-sizeBC/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-3')
    ## FACE-BC-4
    findat=(spanX-sizeBC/2, spanY-sizeBC/2, 0.0)
    modpa.Set(faces=modpa.faces.findAt((findat, ), ), name='FACE-BC-4')
del findat
## BC-CELLS (for frozen area)
if sizeBC!=0 and nely!=1:
    findat1=(sizeBC/2, sizeBC/2, float(height-minSizeTop)/2)
    findat2=(spanX-sizeBC/2, sizeBC/2, float(height-minSizeTop)/2)
    findat3=(sizeBC/2, spanY-sizeBC/2, float(height-minSizeTop)/2)
    findat4=(spanX-sizeBC/2, spanY-sizeBC/2, float(height-minSizeTop)/2)
    modpa.Set(cells=modpa.cells.findAt((findat1, ), (findat2, ), (findat3, ), (findat4, ), ), 
        name='BC-CELLS')
    del findat1, findat2, findat3, findat4

# Divide each cell of the boundary condition in four cells
if sizeBC!=0 and nely==1:
    ## X-DIRECTION
    for x in [sizeBC/2, spanX-sizeBC/2]:
        for y in [meshSize/2]:
            ## Define cell, edge and point
            fCell=(x-0.1, y, (height-minSizeTop)/2)
            fPoint=(x, y-meshSize/2, 0.0)
            fEdge=(x, y-meshSize/2, 0.0)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge    
if sizeBC!=0 and nely!=1:
    ## X-DIRECTION
    for x in [sizeBC/2, spanX-sizeBC/2]:
        for y in [sizeBC/2, spanY-sizeBC/2]:
            ## Define cell, edge and point
            fCell=(x-0.1, y, (height-minSizeTop)/2)
            fPoint=(x, y-sizeBC/2, 0.0)
            fEdge=(x, y-sizeBC/2, 0.0)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge
    ## Y-DIRECTION
    for x in [sizeBC/2, sizeBC, spanX-sizeBC/2, spanX]:
        for y in [sizeBC/2, spanY-sizeBC/2]:
            ## Define cell, edge and point
            fCell=(x-0.1, y, (height-minSizeTop)/2)
            fPoint=(x, y, 0.0)
            fEdge=(x, y+0.1, 0.0)
            ## Apply functions
            partitionEdge(fPoint)
            partitionCell(fCell, fEdge, fPoint)
            del fCell, fPoint, fEdge


########################################################################################################
### SETTING CONSTRAINTS FOR BOUNDARY CONDITIONS                                                      ###
########################################################################################################

# Create reference points and sets
if nely==1 and sizeBC!=0:
    ## RP-1
    findat=(sizeBC/2, meshSize/2, 0.0)
    modra.ReferencePoint(point=modrain.InterestingPoint(modrain.edges.findAt(findat, ), MIDDLE))
    listRP=listRP=modra.referencePoints.keys()
    numRP1=listRP[0]
    modra.Set(name='RP-1', referencePoints=(modra.referencePoints[numRP1], ))
    ## RP-2
    findat=(spanX-sizeBC/2, meshSize/2, 0.0)
    modra.ReferencePoint(point=modrain.InterestingPoint(modrain.edges.findAt(findat, ), MIDDLE))
    listRP=listRP=modra.referencePoints.keys()
    mySet=set(listRP)-set([numRP1])
    numRP2=list(mySet)[0]
    modra.Set(name='RP-2', referencePoints=(modra.referencePoints[numRP2], ))
    del mySet
if nely!=1 and sizeBC!=0:
    ## RP-1
    findat=(sizeBC/2, sizeBC/2, 0.0)
    modra.ReferencePoint(point=modrain.vertices.findAt(findat, ))
    listRP=listRP=modra.referencePoints.keys()
    numRP1=listRP[0]
    modra.Set(name='RP-1', referencePoints=(modra.referencePoints[numRP1], ))
    ## RP-2
    findat=(spanX-sizeBC/2, sizeBC/2, 0.0)
    modra.ReferencePoint(point=modrain.vertices.findAt(findat, ))
    listRP=listRP=modra.referencePoints.keys()
    mySet=set(listRP)-set([numRP1])
    numRP2=list(mySet)[0]
    modra.Set(name='RP-2', referencePoints=(modra.referencePoints[numRP2], ))
    ## RP-3
    findat=(sizeBC/2, spanY-sizeBC/2, 0.0)
    modra.ReferencePoint(point=modrain.vertices.findAt(findat, ))
    listRP=listRP=modra.referencePoints.keys()
    mySet=set(listRP)-set([numRP1, numRP2])
    numRP3=list(mySet)[0]
    modra.Set(name='RP-3', referencePoints=(modra.referencePoints[numRP3], ))
    ## RP-4
    findat=(spanX-sizeBC/2, spanY-sizeBC/2, 0.0)
    modra.ReferencePoint(point=modrain.vertices.findAt(findat, ))
    listRP=listRP=modra.referencePoints.keys()
    mySet=set(listRP)-set([numRP1, numRP2, numRP3])
    numRP4=list(mySet)[0]
    modra.Set(name='RP-4', referencePoints=(modra.referencePoints[numRP4], ))
    del mySet

# Create sets for boundary conditions: if BC is a node and not a surface
if sizeBC==0:
    ## NODE-BC-1
    findat=(0.0, 0.0, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-BC-1')
    ## NODE-BC-2
    findat=(spanX, 0.0, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-BC-2')
    ## NODE-BC-3
    findat=(0.0, spanY, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-BC-3')
    ## NODE-BC-4
    findat=(spanX, spanY, 0.0)
    modpa.Set(vertices=modpa.vertices.findAt((findat, )), name='NODE-BC-4')

# Create constraints (rigid body): if BC is a surface
if sizeBC!=0:
    for i in range(1, len(listRP)+1):
        mod.RigidBody(name='CONSTRAINT-BC-{}'.format(i), pinRegion=modrain.sets['FACE-BC-{}'.format(i)], 
            refPointRegion=modra.sets['RP-{}'.format(i)])

        
########################################################################################################
### SET BOUNDARY CONDITIONS                                                                          ###
########################################################################################################

# Clean script
moddi=mod.DisplacementBC

# Set boundary conditions
if nely==1 and sizeBC!=0:
    ## RP-1
    moddi(createStepName='Initial', name='BC-RP-1', region=modra.sets['RP-1'], u1=SET, u2=SET, u3=SET)
    ## RP-2
    moddi(createStepName='Initial', name='BC-RP-2', region=modra.sets['RP-2'], u2=SET, u3=SET)
if nely!=1 and sizeBC!=0:
    ## RP-1
    moddi(createStepName='Initial', name='BC-RP-1', region=modra.sets['RP-1'], u2=SET, u3=SET)
    ## RP-2
    moddi(createStepName='Initial', name='BC-RP-2', region=modra.sets['RP-2'], u2=SET, u3=SET)
    ## RP-3
    moddi(createStepName='Initial', name='BC-RP-3', region=modra.sets['RP-3'], u1=SET, u3=SET)
    ## RP-4
    moddi(createStepName='Initial', name='BC-RP-4', region=modra.sets['RP-4'], u3=SET)
if sizeBC==0:
    ## BC-NODE-1
    moddi(createStepName='Initial', name='BC-NODE-1', region=modrain.sets['NODE-BC-1'], u2=SET, u3=SET)
    ## BC-NODE-2
    moddi(createStepName='Initial', name='BC-NODE-2', region=modrain.sets['NODE-BC-2'], u2=SET, u3=SET)
    ## BC-NODE-3
    moddi(createStepName='Initial', name='BC-NODE-3', region=modrain.sets['NODE-BC-3'], u1=SET, u3=SET)
    ## BC-NODE-4
    moddi(createStepName='Initial', name='BC-NODE-4', region=modrain.sets['NODE-BC-4'], u3=SET)

# Boundary condtiion in y-direction for 2D problem
if fixy==True:
    ## BC-Y-1
    moddi(createStepName='Initial', name='FACE-Y-1', region=modrain.sets['FACE-Y-1'], u2=SET)
    ## BC-Y-2
    moddi(createStepName='Initial', name='FACE-Y-2', region=modrain.sets['FACE-Y-2'], u2=SET)


########################################################################################################
### CREATING STEP, LOADS AND MESH                                                                    ###
########################################################################################################

# Create step
mod.StaticStep(name='STEP-1', previous='Initial')

# Set field output
del mod.fieldOutputRequests['F-Output-1']
mod.FieldOutputRequest(createStepName='STEP-1', name='F-Output-1', 
    variables=('S', 'MISES', 'E', 'EE', 'NE', 'U', 'RF', 'P', 'COORD', 'ENER', 'ELEN', 'ELEDEN'))

# Set element type
if elType=='C3D8':
    modpa.setElementType(elemTypes=(ElemType(elemCode=C3D8, distortionControl=DEFAULT), 
        ElemType(elemCode=C3D6, distortionControl=DEFAULT), ElemType(elemCode=C3D4, 
        distortionControl=DEFAULT)), regions=modpa.sets['PART-1'])
elif elType=='C3D8R':
    modpa.setElementType(elemTypes=(ElemType(elemCode=C3D8R, distortionControl=DEFAULT), 
        ElemType(elemCode=C3D6, distortionControl=DEFAULT), ElemType(elemCode=C3D4, 
        distortionControl=DEFAULT)), regions=modpa.sets['PART-1'])
elif elType=='C3D20':
    modpa.setElementType(elemTypes=(ElemType(elemCode=C3D20), ElemType(elemCode=C3D15), 
        ElemType(elemCode=C3D10)), regions=modpa.sets['PART-1'])
elif elType=='C3D20R':
    modpa.setElementType(elemTypes=(ElemType(elemCode=C3D20R), ElemType(elemCode=C3D15), 
        ElemType(elemCode=C3D10)), regions=modpa.sets['PART-1'])
else:
    sys.exit()

# Seed
modpa.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=meshSize)

# Generate mesh
modpa.generateMesh()

# Create job (first check the model before applying topology optimization)
mdb.Job(model=modelName, name=jobName, numCpus=numCores, numDomains=numCores, memory=perMem, 
    memoryUnits=PERCENTAGE)

# Submit job
if submitJob==True:
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()


########################################################################################################
### PREPARATION FOR TOPOLOGY OPTIMIZATION                                                            ###
########################################################################################################

# Settig coordinate system for all models (for planar symmetry)
modra.DatumCsysByThreePoints(coordSysType=CARTESIAN, 
    line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0), name='FLOOR-CENTER', 
    origin=(spanX/2, spanY/2, height/2))
listCsys=modra.datums.keys()
cSysNum=listCsys[-1]


########################################################################################################
### COPY MODELS (CREATE THE DIFFERENT LOAD CASES)                                                    ###
########################################################################################################

if nLoadCases>1:
    for num in range(2, nLoadCases+1):
        mdb.Model(name='LC-{}'.format(num), objectToCopy=mdb.models[modelName])


########################################################################################################
### SET LOADING FOR THE DIFFERENT LOAD CASES                                                         ###
########################################################################################################

# Set pressure load on tiles in top-layer
for lc in range(1, nLoadCases+1):
    mod=mdb.models['LC-{}'.format(lc)]
    modra=mod.rootAssembly
    for tile in range(1, (nTilesX*nTilesY)+1):
        if globals()["LC{}_T{}".format(lc, tile)]!=0:
            mod.Pressure(createStepName='STEP-1', name='LOAD-TILE-{}'.format(tile), 
                magnitude=globals()["LC{}_T{}".format(lc, tile)], distributionType=TOTAL_FORCE, 
                region=modrain.surfaces['LOAD-TILE-{}'.format(tile)])


########################################################################################################
### CREATE TOPOLOGY OPTIMIZATION TASK                                                                ###
########################################################################################################

# Set back values to make optiimization task
mod=mdb.models[modelName]
modra=mod.rootAssembly

# Task
mod.TopologyTask(freezeBoundaryConditionRegions=ON, freezeLoadRegions=ON, 
    materialInterpolationTechnique=SIMP, materialInterpolationPenalty=interpolationPenalty, 
    name=taskName, region=modra.sets['DESIGN-SPACE'], densityMoveLimit=maxDeltaDens, 
    initialDensity=initialDens, maxDensity=maxDens, minDensity=minDens, 
    objectiveFunctionDeltaStopCriteria=conObjective, elementDensityDeltaStopCriteria=conDensity,
    numFulfilledStopCriteria=2) 

# Cleaning script
modopt=mod.optimizationTasks[taskName]

# Design response: strain energy
## Make list with load cases (models) and steps
LC1=(ALL, '', ALL, ALL, '') # Optimization task created in model LC-1
tupleLC=(LC1, )
if nLoadCases>1:
    for num in range(2, nLoadCases+1):
        LCx=(ALL, '', ALL, ALL, 'LC-{}'.format(num))
        tupleLC=tupleLC+(LCx, )
        del LCx
## Define operation
if operation=='SUM':
    ope=SUM
if operation=='MAX':
    ope=MAXIMUM
if operation=='MIN':
    ope=MINIMUM
## Set design respons
modopt.SingleTermDesignResponse(identifier='STRAIN_ENERGY', name='D-RESP-1-STRAIN', operation=SUM, 
    stepOptions=(tupleLC), region=MODEL, stepOperation=ope)

# Design response
## Volume
modopt.SingleTermDesignResponse(drivingRegion=None, operation=SUM, identifier='VOLUME', 
    name='D-RESP-2-VOLUME', stepOptions=(), region=MODEL)
## Von Mises stress
if maxStress!=None:
    modopt.SingleTermDesignResponse(
        name='D-RESP-3-STRESS', region=MODEL, identifier='SIG_TOPO_MISES', 
        drivingRegion=None, operation=MAXIMUM, stepOptions=())

# Objective: strain
modopt.ObjectiveFunction(name='OBJECTIVE-1-STRAIN', 
    objectives= ((OFF, 'D-RESP-1-STRAIN', 1.0, 0.0, ''), ))

# Constraints
## Volume
modopt.OptimizationConstraint(designResponse='D-RESP-2-VOLUME', name='CONSTRAINT-1-VOLUME', 
    restrictionValue=volumeFraction, restrictionMethod=RELATIVE_LESS_THAN_EQUAL)
## Von Mises stress
if maxStress!=None:
    modopt.OptimizationConstraint(name='CONSTRAINT-2-STRESS', designResponse='D-RESP-3-STRESS', 
        restrictionMethod=ABSOLUTE_LESS_THAN_EQUAL, restrictionValue=maxStress)

# Geometry restrictions
## Frozen area (at BC)
if sizeBC!=0 and nely!=1 and freezeCellsBC==True:
    modopt.FrozenArea(name='FROZEN-AREA-BC-CELLS', region=modrain.sets['BC-CELLS'])
## Frozen area (top layer)
if freezeTop==True:
    modopt.FrozenArea(name='FROZEN-AREA-TOP-LAYER', region=modrain.sets['TOP-LAYER'])
## Minimum member size
minSize=2*rmin
if minSize>0:
    mod.optimizationTasks[taskName].TopologyMemberSize(minThickness=minSize, name='MIN-THICKNESS',  
        sizeRestriction=MINIMUM, region=modra.sets['DESIGN-SPACE'])
## Planar symmetry
if planarSym1==True:
    modopt.TopologyPlanarSymmetry(axis=AXIS_1, name='PLANAR-SYM1', region=modra.sets['DESIGN-SPACE'], 
        csys=modra.datums[cSysNum])
if planarSym2==True:
    modopt.TopologyPlanarSymmetry(axis=AXIS_2, name='PLANAR-SYM2', region=modra.sets['DESIGN-SPACE'], 
        csys=modra.datums[cSysNum])
## Demold control
if demold==True:
    p1=(0.0, 0.0, 0.0)
    p2=(0.0, 0.0, height-minSizeTop)
    modopt.TopologyDemoldControl(collisionCheckRegion=DEMOLD_REGION, csys=None, draftAngle=demoldAngle, 
        name='DEMOLD-CONTROL', pointRegion=modra.sets['NODE-DEMOLD'], technique=POINT, 
        pullDirection=(modrain.vertices.findAt(p1, ), modrain.vertices.findAt(p2, )), 
        region=modra.sets['DESIGN-SPACE'])
    del p1, p2

# Topology process
mdb.OptimizationProcess(dataSaveFrequency=OPT_DATASAVE_SPECIFY_CYCLE, saveEvery=None, saveFirst=False, 
    saveInitial=True, maxDesignCycle=mDC, model=modelName, name=processName, odbMergeFrequency=2, 
    prototypeJob=optJobName, task=taskName) 
mdb.optimizationProcesses[processName].Job(model=modelName, name=optJobName, numCpus=numCores, 
    numDomains=numCores, memory=perMem, memoryUnits=PERCENTAGE)  

# Optimization process
if submitOpt==True:
    ## Submit optimization process
    mdb.optimizationProcesses[processName].submit()
    ## Wait for completion
    mdb.optimizationProcesses[processName].waitForCompletion()
    ## Combine results
    mdb.CombineOptResults(analysisFieldVariables=ALL, models=ALL, optIter=ALL, steps=('STEP-1', ), 
        optResultLocation=pathProcess)


########################################################################################################
### END OF SCRIPT                                                                                    ###
########################################################################################################

# End
print('END OF SCRIPT')