# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-3 replay file
# Internal Version: 2013_10_09-13.29.32 126623
# Run by johan on Tue Mar  6 17:02:52 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=719.666625976562, 
    height=237.175018310547)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=5.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7673, 
    farPlane=5.66079, width=10.9823, height=6.01385, cameraPosition=(0.259537, 
    -0.235822, 4.71405), cameraTarget=(0.259537, -0.235822, 0))
s.rectangle(point1=(0.0, 0.0), point2=(5.0, 5.0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=0.05)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Part-1']
p.seedPart(size=0.25, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-1']
s1 = p.features['Solid extrude-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s1)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Solid extrude-1'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=13.1778, 
    farPlane=15.0565, width=9.88728, height=5.80803, cameraPosition=(2.76871, 
    2.50559, 14.1421), cameraTarget=(2.76871, 2.50559, 0))
s2.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
del p.features['Solid extrude-1']
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=5000.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3823.34, 
    farPlane=5604.75, width=9631.59, height=5657.83, cameraPosition=(-434.005, 
    -46.3025, 4714.05), cameraTarget=(-434.005, -46.3025, 0))
s.rectangle(point1=(-2500.0, -1000.0), point2=(2500.0, 4000.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3429.41, 
    farPlane=5998.68, width=13341.3, height=7836.99, cameraPosition=(-98.9064, 
    367.38, 4714.05), cameraTarget=(-98.9064, 367.38, 0))
p = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-2']
p.BaseSolidExtrude(sketch=s, depth=50.0)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Model-1'].parts['Part-1']
p = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-2']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-2']
s1 = p.features['Solid extrude-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s1)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Solid extrude-1'], filter=COPLANAR_EDGES)
s2.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].setValues(displayedObject=a)
#* FeatureError: Regeneration failed
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a = mdb.models['Model-1'].rootAssembly
del a.features['Part-1-1']
a1 = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-2']
a1.Instance(name='Part-2-1', part=p, dependent=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p1 = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['Part-2']
p.seedPart(size=50.0, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-2']
p.seedPart(size=100.0, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-2']
p.generateMesh()
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['Part-2']
c = p.cells
pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-2']
e = p.edges
pickedEndEdges = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=pickedEndEdges, minSize=5.0, 
    maxSize=300.0, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-2']
e = p.edges
pickedCenterEdges = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=pickedCenterEdges, minSize=5.0, 
    maxSize=300.0, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-2']
e = p.edges
pickedEndEdges = e.getSequenceFromMask(mask=('[#200 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=pickedEndEdges, minSize=5.0, 
    maxSize=300.0, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-2']
e = p.edges
pickedCenterEdges = e.getSequenceFromMask(mask=('[#200 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=pickedCenterEdges, minSize=5.0, 
    maxSize=300.0, constraint=FINER)
session.viewports['Viewport: 1'].view.setValues(nearPlane=10985.3, 
    farPlane=17749.6, width=8282.55, height=4878.16, cameraPosition=(-8099.3, 
    8086.51, 9896.65), cameraUpVector=(0.43971, 0.688775, -0.576406), 
    cameraTarget=(128.382, 1497.97, -101.348))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10958.6, 
    farPlane=17741.8, width=8262.44, height=4866.31, cameraPosition=(-10745.2, 
    5690.48, 8564.01), cameraUpVector=(0.375147, 0.798182, -0.471349), 
    cameraTarget=(157.954, 1524.75, -86.4538))
p = mdb.models['Model-1'].parts['Part-2']
e = p.edges
pickedCenterEdges = e.getSequenceFromMask(mask=('[#5 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=pickedCenterEdges, minSize=5.0, 
    maxSize=300.0, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-2']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=11106.8, 
    farPlane=17593.5, width=6845.23, height=4031.62, viewOffsetX=-1970.59, 
    viewOffsetY=-250.598)
session.viewports['Viewport: 1'].view.setValues(nearPlane=11068.3, 
    farPlane=17578, width=6821.45, height=4017.62, cameraPosition=(-10745.2, 
    5690.48, 8518.62), cameraTarget=(157.954, 1524.75, -131.844), 
    viewOffsetX=-1963.74, viewOffsetY=-249.727)
session.viewports['Viewport: 1'].view.setValues(nearPlane=10842.3, 
    farPlane=15701.3, width=6682.2, height=3935.6, cameraPosition=(-5368.17, 
    5690.48, 11424.4), cameraUpVector=(0.133772, 0.798182, -0.587376), 
    cameraTarget=(708.047, 1524.75, -1097.1), viewOffsetX=-1923.65, 
    viewOffsetY=-244.629)
session.viewports['Viewport: 1'].view.setValues(nearPlane=10750, 
    farPlane=15793.5, width=7946.46, height=4680.22, viewOffsetX=-2386.43, 
    viewOffsetY=-431.65)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(-46.1045, 
    16236.7, -37.2188), cameraUpVector=(0, 0, 1))
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(14570.8, 
    1619.77, -37.2188), cameraUpVector=(0, 1, 0))
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(-46.1045, 
    1619.77, 14579.7))
session.viewports['Viewport: 1'].view.setValues(nearPlane=13702.2, 
    farPlane=15407.2, width=8913.89, height=5250, viewOffsetX=-126.863, 
    viewOffsetY=-37.0503)
mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.Job(name='3350elm-dense-center', model='Model-1', description='', 
    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
mdb.jobs['3350elm-dense-center'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "3350elm-dense-center.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(128.382, 
    1497.96, 14426.7), cameraUpVector=(0, 1, 0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=13409.3, 
    farPlane=15394.1, width=10136.7, height=5954.56, cameraPosition=(14.1759, 
    1497.96, 14426.7), cameraTarget=(14.1763, 1497.96, -101.348))
session.viewports['Viewport: 1'].view.setValues(nearPlane=13419.4, 
    farPlane=15384, width=10144.3, height=5959.04, cameraPosition=(0.880002, 
    1497.96, 14426.7), cameraTarget=(0.880402, 1497.96, -101.348))
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(14528.9, 
    1497.96, -101.348))
session.viewports['Viewport: 1'].view.setValues(nearPlane=11239.1, 
    farPlane=17818.7, width=8496.17, height=4990.86, cameraPosition=(14528.9, 
    1536.41, -101.348), cameraTarget=(0.880402, 1536.41, -101.348))
session.viewports['Viewport: 1'].view.setValues(nearPlane=12001.1, 
    farPlane=17251.2, width=9072.21, height=5329.24, cameraPosition=(-6137.5, 
    -2802.96, -12534.8), cameraUpVector=(-0.248192, 0.94618, -0.20771), 
    cameraTarget=(-0.371918, 1536.15, -102.101))
session.viewports['Viewport: 1'].view.setValues(nearPlane=11382.6, 
    farPlane=17838.1, width=8604.66, height=5054.59, cameraPosition=(-10602, 
    -1571.43, -9547.59), cameraUpVector=(-0.129255, 0.975718, -0.176825), 
    cameraTarget=(-30.3159, 1544.41, -82.0654))
session.viewports['Viewport: 1'].view.setValues(nearPlane=11423.6, 
    farPlane=17797.1, width=8635.66, height=5072.8, cameraPosition=(-10602, 
    -1571.43, -9547.59), cameraUpVector=(-0.139298, 0.976273, -0.165791), 
    cameraTarget=(-30.3159, 1544.41, -82.0654))
