function xieta = inverseIsoparametricTransformation(x,Xe,nodesCoord,degree)

global degreeGlobal XeGlobal nodesCoordGlobal xgoalGlobal
degreeGlobal = degree; XeGlobal=Xe; nodesCoordGlobal = nodesCoord;

tol = 1.e-10;

%First trial assuming straight-sided element
xieta0 = inverseLinearTransformation(x,Xe); 
x0 = isoparametricTransformationHighOrder(xieta0,Xe,degree,nodesCoord);

if norm(x-x0)<tol*norm(x)+1.e-14
    xieta = xieta0;
else
    xgoalGlobal=x;
    xieta = fzero(@fnonlinearInverseIsoTrans,xieta0);
end

    

    

