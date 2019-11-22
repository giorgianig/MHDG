function f=fnonlinearInverseIsoTrans(xieta)

global degreeGlobal XeGlobal nodesCoordGlobal xgoalGlobal

x=isoparametricTransformationHighOrder(xieta,XeGlobal,degreeGlobal,nodesCoordGlobal);
f = x-xgoalGlobal;

