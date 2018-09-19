########################################################################################################
### 2D BEAM: BENCHMARK FOR TOPOLOGY OPTIMIZATION                                                     ###
### Simply supported beam                                                                            ###
### Arjen Deetman                                                                                    ###
### version 17-09-2018                                                                               ###
########################################################################################################

########################################################################################################
### INPUT PARAMETERS                                                                                 ###
########################################################################################################

# File and job names
jobName='MODEL-CHECK' ## Job name to check the model
taskName='OPTIMIZATION-TASK' ## Task name for topology optimization task
processName='TEST-OPT-2D-MODEL' ## Process name for topology optimization
optJobName='{}-JOB'.format(processName) ## Job name for topology optimization

# Mesh szize
meshSize=5.0 ## Mesh size in mm

# Geometry
witdth=float(1200) ## Width in mm
height=float(200) ## Height in mm
LS=float(40) ## Length of support in mm 

# Geometry restrictions
rmin=float(1.5*meshSize) ## Filter radius
minSizeTop=float(meshSize) ## Minimum thickness of the top floor in mm
planarSym=True ## Planar symmetry (for axis 1)
demold=False ## Demold control

# Material properties
E=float(10000) ## E-modulus in N/mm2
p=float(0.2) ## Poison ratio 

# Optimization parameters
VF=float(0.5) ## Volume fraction (for design area)
operation='SUM' ## MAX, MIN or SUM: Operation between steps and load cases (only for multiple lc)

# Optimization task settings
## Material interpolation penalty
interpolationPenalty=3.0 ## Material interpolation penalty (default=3.0)
## Density
minDens=0.001 ## Minimum density (default=0.001)
maxDens=1.0 ## Maximum density (default=0.001)
maxDeltaDens=0.20 ## Maximum change per design cycle (default=0.25)
initialDens=VF ## Initial density (default=DEFAULT)
## Convergence
### Stops when both criteria are reached. 
conObjective=float(0.001) ## Objective function delta criterion (default=0.001)
conDensity=float(0.005) ## Element density delta criterion (default=0.005)
## Cycles
mDC=int(1000) ## Maximum number of cycles

# Computer settings
numCores=int(2) # Number of cores that will be uese of the processor
perMem=int(90) # Percentage of the memory that will be uses of the computer in %

# Submit
submitJob=False ## Normal job to check the model (only for LC-1)
submitOpt=False ## Submitting the optimization task


########################################################################################################
### DEFINE LOADCASES AND LOADS                                                                       ###
########################################################################################################

# Number of load cases
nLoadCases=int(1)

# Number of tiles in top layer
nTiles=int(60) # Number of edges to apply a load (top will be divided in n)

# Settings for all load cases
tilePressure=float(1) ## Pressure

# Load
## Set load for all edges in top layer
for lc in range(1, nLoadCases+1):
    for tile in range(1, nTiles+1):
        globals()["LC{}_T{}".format(lc, tile)]=0.0
## Load case 1: line load (pressure on all create tiles)
for tile in range(1, nTiles+1):
    globals()["LC{}_T{}".format(lc, tile)]=tilePressure


########################################################################################################
### CREATE ABAQUS MODEL                                                                              ###
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

# Clean script
mod=mdb.models[modelName]

# x,y
x=float(witdth*0.5)
y=float(height*0.5)

# Create part 
mod.ConstrainedSketch(name='__profile__', sheetSize=200.0)
## Cleaning script
modske=mod.sketches['__profile__']
## Vertical left
### Line 
p1,p2=(-x, y),(-x, -y+LS)
modske.Line(point1=p1, point2=p2)
### Line
p1,p2=p2,(-x,-y)
modske.Line(point1=p1, point2=p2) 
## Horizontal bottom 
### Line
p1,p2=p2,(-x+LS, -y)
modske.Line(point1=p1, point2=p2) 
### Line
p1,p2=p2,(0.0, -y)
modske.Line(point1=p1, point2=p2) 
### Line
p1,p2=p2,(x-LS, -y)
modske.Line(point1=p1, point2=p2) 
### Line
p1,p2=p2,(x, -y)
modske.Line(point1=p1, point2=p2) 
## Vertical right 
### Line
p1,p2=p2,(x, -y+LS)
modske.Line(point1=p1, point2=p2) 
### Line
p1,p2=p2,(x, y)
modske.Line(point1=p1, point2=p2)
### Parameters for load edges
LL=float(2*x)/nTiles 
### Horizontal upper (divided in n lines)
for a in range(0, nTiles+1):
    p1,p2=p2,(x-a*LL, y)
    modske.Line(point1=p1, point2=p2) 
## Create part 
mod.Part(dimensionality=TWO_D_PLANAR, name='PART-1', type=DEFORMABLE_BODY) 
mod.parts['PART-1'].BaseShell(sketch=modske) 
del modske, mod.sketches['__profile__'], p1, p2

# Clean script
modpa=mod.parts['PART-1']

# Create set of the whole model: PART-1
modpa.Set(faces=modpa.faces.getSequenceFromMask(('[#1 ]', ), ), name='PART-1') 

# Create material
mod.Material(name='STEEL') 
mod.materials['STEEL'].Elastic(table=((E, p), )) 

# Create and section
mod.HomogeneousSolidSection(material='STEEL', name='SECTION-1', thickness=None) 

# Section Assignment 
modpa.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
    region=modpa.sets['PART-1'], sectionName='SECTION-1', thicknessAssignment=FROM_SECTION) 
 
# Cleaning script
modra=mod.rootAssembly

# Assembly 
modra.DatumCsysByDefault(CARTESIAN) 
modra.Instance(dependent=OFF, name='PART-1-1', part=modpa) 

# Clean script
modrain=modra.instances['PART-1-1']

# Creating partition: split in bottom and top
## Top is minimum thickness of topfloor
## Bottom is design area for topology optimization
if minSizeTop>0:
    mod.ConstrainedSketch(gridSpacing=100.0, name='__profile__', 
        sheetSize=4000.0, transform=modpa.MakeSketchTransform(
        sketchPlane=modpa.faces.findAt((.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    modske=mod.sketches['__profile__']
    modpa.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=modske)
    mod.sketches['__profile__'].Line(point1=(-x, y-minSizeTop), point2=(x, y-minSizeTop))
    modpa.PartitionFaceBySketch(faces=modpa.faces.findAt(((0.0, 0.0, 0.0), )), sketch=modske)
    del modske, mod.sketches['__profile__']
    modra.regenerate()

# Create sets
## Top layer (minimum thickness of top floor)
if minSizeTop>0:
    findat=(0.0, y-float(minSizeTop)/2, 0.0)
    modra.Set(faces=modrain.faces.findAt((findat, )), name='TOP-LAYER')
## Bottom layer
findat=(0.0, y-minSizeTop-float(y*2-minSizeTop)/2, 0.0)
modra.Set(faces=modrain.faces.findAt((findat, )), name='BOTTOM-LAYER')
## Design area
findat1=(0.0, y-float(minSizeTop)/2, 0.0)
findat2=(0.0, y-minSizeTop-float(y*2-minSizeTop)/2, 0.0)
modra.Set(faces=modrain.faces.findAt((findat1, ), (findat2, )), name='DESIGN-AREA')
del findat, findat1, findat2
## Edge for applying load
for a in range(0, nTiles):
    findat=(-x+0.5*LL+a*LL, y, 0.0)
    modra.Surface(name='LOAD-TILE-{}'.format(a+1), side1Edges=modrain.edges.findAt((findat, ))) 
## Edge for applying support 1
findat=(-x+0.5*LS, -y, 0.0)
modra.Set(name='EDGE-BC-1', edges=modrain.edges.findAt((findat, )))
## Edge for applying support 2
findat=(x-0.5*LS, -y, 0.0)
modra.Set(name='EDGE-BC-2', edges=modrain.edges.findAt((findat, ))) 
## Node for demold control
findat=(-x, y-minSizeTop, 0.0)
modra.Set(name='NODE-DEMOLD', vertices=modrain.vertices.findAt((findat, )))
del findat

# Create step 
mod.StaticStep(name='STEP-1', previous='Initial') 

# Set field output
## Delete standard fiel output
del mod.fieldOutputRequests['F-Output-1']
## Create new field output
mod.FieldOutputRequest(createStepName='STEP-1', name='F-Output-1', 
    variables=('S', 'MISES', 'E', 'EE', 'NE', 'U', 'RF', 'P', 'COORD', 'ENER', 'ELEN', 'ELEDEN'))

# Create reference point for boundary conditions
## RF-1
findat=(-x+0.5*LS, -y, 0.0)
modra.ReferencePoint(point=modrain.InterestingPoint(modrain.edges.findAt(findat, ), MIDDLE))
listRP=listRP=modra.referencePoints.keys()
numRP1=listRP[0]
modra.Set(name='RP-1', referencePoints=(modra.referencePoints[numRP1], ))
## RF-2
findat=(x-0.5*LS, -y, 0.0)
modra.ReferencePoint(point=modrain.InterestingPoint(modrain.edges.findAt(findat, ), MIDDLE))
listRP=listRP=modra.referencePoints.keys()
mySet=set(listRP)-set([numRP1])
numRP2=list(mySet)[0]
modra.Set(name='RP-2', referencePoints=(modra.referencePoints[numRP2], ))
del mySet

# Create rigid body constraints
for i in [1,2]:
    mod.RigidBody(name='CONSTRAINT-BC-{}'.format(i), 
        pinRegion=modra.sets['EDGE-BC-{}'.format(i)], 
        refPointRegion=modra.sets['RP-{}'.format(i)])

# Set boundary conditions
## RP-1
mod.DisplacementBC(createStepName='Initial', name='BC-RP-1', region=modra.sets['RP-1'], u1=SET, u2=SET)
## RP-2
mod.DisplacementBC(createStepName='Initial', name='BC-RP-2', region=modra.sets['RP-2'], u2=SET)

# Set mesh size 
modra.seedPartInstance(deviationFactor=0.1,  minSizeFactor=0.1, regions=(modrain, ), size=meshSize) 

# Generate mesh
modra.generateMesh(regions=(modrain, )) 

# Create job (first check the model before applying topology optimization)
mdb.Job(model=modelName, name=jobName, numCpus=numCores, numDomains=numCores,  memory=perMem, 
    memoryUnits=PERCENTAGE)
if submitJob==True:
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()


########################################################################################################
### COPY MODELS (CREATE THE DIFFERENT LOAD CASES)                                                    ###
########################################################################################################

if nLoadCases>1:
    for num in range(2, nLoadCases+1):
        mdb.Model(name='LC'.format(num), objectToCopy=mdb.models[modelName])


########################################################################################################
### SET LOADING FOR THE DIFFERENT LOAD CASES                                                         ###
########################################################################################################

# Set load 
for lc in range(1, nLoadCases+1):
    mod=mdb.models['LC-{}'.format(lc)]
    modra=mod.rootAssembly
    for tile in range(1, (nTiles+1)):
        if globals()["LC{}_T{}".format(lc, tile)]!=0:
            mod.Pressure(createStepName='STEP-1', region=modra.surfaces['LOAD-TILE-{}'.format(tile)], 
                name='LOAD-TILE-{}'.format(tile), magnitude=globals()["LC{}_T{}".format(lc, tile)]) 


########################################################################################################
### TOPOLOGY OPTIMIZATION                                                                            ###
########################################################################################################

# Set back values to make optiimization task
mod=mdb.models[modelName]
modra=mod.rootAssembly

# Task
mod.TopologyTask(freezeBoundaryConditionRegions=ON, freezeLoadRegions=ON, 
    materialInterpolationTechnique=SIMP, materialInterpolationPenalty=interpolationPenalty, 
    name=taskName, region=modra.sets['DESIGN-AREA'], densityMoveLimit=maxDeltaDens, 
    initialDensity=initialDens, maxDensity=maxDens, minDensity=minDens, 
    objectiveFunctionDeltaStopCriteria=conObjective, elementDensityDeltaStopCriteria=conDensity) 

# Cleaning script
modopt=mod.optimizationTasks[taskName]

# Design response: strain energy
## Make list with load cases (models) and steps
LC1=(ALL, '', ALL, ALL, '') # Optimization task created in LC-1
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
modopt.SingleTermDesignResponse(drivingRegion=None, identifier='VOLUME', name='D-RESP-2-VOLUME', 
    operation=SUM, region=MODEL, stepOptions=())

# Objective: strain
modopt.ObjectiveFunction(name='OBJECTIVE-1-STRAIN', 
    objectives=((OFF, 'D-RESP-1-STRAIN', 1.0, 0.0, ''), ))

# Constraints
## Volume
modopt.OptimizationConstraint(designResponse='D-RESP-2-VOLUME', name='CONSTRAINT-1-VOLUME', 
    restrictionValue=VF, restrictionMethod=RELATIVE_LESS_THAN_EQUAL)

# Geometry restrictions
## Frozen area (top layer)
if minSizeTop>0:
    modopt.FrozenArea(name='FROZEN-AREA-TOP-LAYER', region=modra.sets['TOP-LAYER'])
## Minimum member size
minSize=2*rmin
if minSize>0:
    mod.optimizationTasks[taskName].TopologyMemberSize(minThickness=minSize, 
        region=modra.sets['DESIGN-AREA'], sizeRestriction=MINIMUM, name='MIN-THICKNESS')
## Planar symmetry
if planarSym==True:
    modopt.TopologyPlanarSymmetry(axis=AXIS_1, name='PLANAR-SYM', region=modra.sets['DESIGN-AREA'])
## Demold control
if demold==True:
    modopt.TopologyDemoldControl(collisionCheckRegion=DEMOLD_REGION, csys=None, draftAngle=0.0, 
        name='DEMOLD-CONTROL', pointRegion=modra.sets['NODE-DEMOLD'],  technique=POINT,
        pullDirection=(modrain.vertices[3], modrain.vertices[1]), region=modra.sets['BOTTOM-LAYER'])

# Topology process 
mdb.OptimizationProcess(dataSaveFrequency=OPT_DATASAVE_SPECIFY_CYCLE, saveEvery=None, saveFirst=False, 
    saveInitial=True, maxDesignCycle=mDC, model=modelName, name=processName, odbMergeFrequency=2, 
    prototypeJob=optJobName, task=taskName) 
mdb.optimizationProcesses[processName].Job(model=modelName, name=optJobName, 
    numCpus=numCores, numDomains=numCores, memory=perMem, memoryUnits=PERCENTAGE) 

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