%load HDGCoarse_soln
load HDGFine_soln

x = [0.1,2.6]; 
x = [8,-2.4];
x = [8,-3];
%x = X(T(103,1),:);

figure(1), clf, plotMesh(X,T); hold on, plot(x(1),x(2),'ro','LineWidth',2)
[elementNum,xieta] = findElementForx(x,X,T,referenceElement.NodesCoord,referenceElement.degree)