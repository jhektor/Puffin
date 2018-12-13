# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-3 replay file
# Internal Version: 2013_10_09-13.29.32 126623
# Run by johan on Thu Mar  8 10:16:46 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=419.099975585938, 
    height=237.175018310547)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
session.viewports['Viewport: 1'].view.setValues(width=1.97746, height=1.16161, 
    viewOffsetX=0.0809054, viewOffsetY=-0.0106618)
openMdb(pathName='/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae')
#: The model database "/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.Set(faces=faces, name='top')
#: The set 'top' has been created (1 face).
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.Set(faces=faces, name='right')
#: The set 'right' has been created (1 face).
mdb.save()
#: The model database has been saved to "/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae".
session.viewports['Viewport: 1'].view.setValues(nearPlane=8802.71, 
    farPlane=16283.6, width=8318.82, height=4209.64, cameraPosition=(5860.61, 
    11220.4, -3306.57), cameraUpVector=(-0.404115, 0.191096, 0.894524), 
    cameraTarget=(110.019, 914.546, 975.436))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9596.38, 
    farPlane=15596.8, width=9068.86, height=4589.19, cameraPosition=(5339.52, 
    4118.96, -9975.06), cameraUpVector=(-0.107419, 0.815493, 0.568711), 
    cameraTarget=(110.49, 920.97, 981.468))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9048.13, 
    farPlane=16027.3, width=8550.75, height=4327, cameraPosition=(-9547.36, 
    -2319.5, -6418.67), cameraUpVector=(0.267115, 0.944387, -0.19179), 
    cameraTarget=(60.744, 899.455, 993.352))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#20 ]', ), )
p.Set(faces=faces, name='back')
#: The set 'back' has been created (1 face).
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#8 ]', ), )
p.Set(faces=faces, name='bottom')
#: The set 'bottom' has been created (1 face).
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.Set(faces=faces, name='left')
#: The set 'left' has been created (1 face).
mdb.save()
#: The model database has been saved to "/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae".
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=8872, 
    farPlane=16039.6, width=8384.3, height=4253.93, cameraPosition=(-9951.11, 
    5131.65, 7249.72), cameraUpVector=(0.223526, 0.634928, -0.739529), 
    cameraTarget=(61.2848, 889.48, 975.053))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10374.8, 
    farPlane=14570.7, width=9804.45, height=4974.47, cameraPosition=(-640.026, 
    2240.67, 13394.6), cameraUpVector=(0.0597296, 0.901484, -0.428672), 
    cameraTarget=(-12.4974, 912.389, 926.361))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10240.4, 
    farPlane=14705, width=9677.42, height=4910.02, cameraPosition=(-640.026, 
    2240.67, 13394.6), cameraUpVector=(-6.69621e-05, 0.902, -0.431736), 
    cameraTarget=(-12.4974, 912.389, 926.361))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9028.03, 
    farPlane=15912.2, width=8531.71, height=4328.72, cameraPosition=(-10590.7, 
    2947.92, -5289.06), cameraUpVector=(0.463814, 0.871037, 0.161773), 
    cameraTarget=(52.7449, 907.752, 1048.86))
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEndEdges = e.getSequenceFromMask(mask=('[#281 ]', ), )
pickedCenterEdges = e.getSequenceFromMask(mask=('[#4 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=pickedEndEdges, 
    centerEdges=pickedCenterEdges, minSize=10.0, maxSize=200.0, 
    constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#fff ]', ), )
p.deleteSeeds(regions=pickedEdges)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9096.03, 
    farPlane=15764.3, width=8595.98, height=4361.33, cameraPosition=(-4281.88, 
    7505.82, -8687.82), cameraUpVector=(0.319526, 0.624649, 0.712542), 
    cameraTarget=(10.0587, 876.913, 1071.86))
session.viewports['Viewport: 1'].view.setValues(nearPlane=8772.95, 
    farPlane=16015.2, width=8290.66, height=4206.42, cameraPosition=(5624.32, 
    8879.58, -6739.27), cameraUpVector=(-0.166874, 0.449733, 0.877436), 
    cameraTarget=(-88.9993, 863.176, 1052.38))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9129.48, 
    farPlane=15692.8, width=8627.59, height=4377.37, cameraPosition=(5847.75, 
    5623.57, -8923.02), cameraUpVector=(-0.409953, 0.744959, 0.526285), 
    cameraTarget=(-91.8906, 905.309, 1080.64))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9299.62, 
    farPlane=15535.1, width=8788.37, height=4458.95, cameraPosition=(5162.31, 
    4743.93, -9654.92), cameraUpVector=(-0.31928, 0.79442, 0.516679), 
    cameraTarget=(-83.9753, 915.467, 1089.09))
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedCenterEdges = e.getSequenceFromMask(mask=('[#285 ]', ), )
p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=pickedCenterEdges, 
    minSize=10.0, maxSize=200.0, constraint=FINER)
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
p.seedPart(size=400.0, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Solid extrude-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Solid extrude-1'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=8208.17, 
    farPlane=12419.2, width=11481.3, height=5809.99, cameraPosition=(-134.344, 
    978.716, 11313.7), cameraTarget=(-134.344, 978.716, 0))
s1.autoDimension(objectList=(g[2], g[3], g[4], g[5]))
#: 2 dimensions added
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
mdb.Job(name='box-4x4x2um', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
mdb.jobs['box-4x4x2um'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "box-4x4x2um.inp".
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Solid extrude-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Solid extrude-1'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=8429.21, 
    farPlane=12198.2, width=9564.55, height=4840.03, cameraPosition=(-333.391, 
    956.169, 11313.7), cameraTarget=(-333.391, 956.169, 0))
s2.autoDimension(objectList=(g[2], g[3], g[4], g[5]))
#: 2 dimensions added
mdb.save()
#: The model database has been saved to "/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae".
s2.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
mdb.save()
#: The model database has been saved to "/home/johan/projects/Puffin/inputfiles/multiSn/abq-mesh.cae".
