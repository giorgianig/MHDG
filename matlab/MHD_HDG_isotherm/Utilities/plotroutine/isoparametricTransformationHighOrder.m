function x=isoparametricTransformationHighOrder(xieta,Xe,degree,nodesCoord) 

shapeFunctions = computeShapeFunctionsAtPoints(degree,nodesCoord,xieta);
N = shapeFunctions(:,:,1)';
x = N*Xe;