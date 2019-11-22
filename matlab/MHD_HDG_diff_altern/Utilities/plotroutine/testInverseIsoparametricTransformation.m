nodesCoord = [   -1    -1
     1    -1
    -1     1
     0    -1
     0     0
    -1     0];
degree = 2;

Xe =[    0    3.0000
         0    2.5000
    0.5000    2.5000
         0    2.7500
    0.2500    2.5000
    0.2500    2.800];
Te = 1:size(Xe,1);

x = [0.1,2.6]; 
figure(1), plotMesh(Xe,Te), hold on, plot(Xe(:,1),Xe(:,2),'k*',x(1),x(2),'ro','LineWidth',2), hold off

xieta0 = inverseLinearTransformation(x,Xe(1:3,:))
figure(2), plotMesh(nodesCoord,Te), hold on, plot(nodesCoord(:,1),nodesCoord(:,2),'k*',xieta0(1),xieta0(2),'ro','LineWidth',2), hold off

x0 = isoparametricTransformationHighOrder(xieta0,Xe,degree,nodesCoord)

xieta = inverseIsoparametricTransformation(x,Xe,nodesCoord,degree)
figure(2), hold on, plot(xieta(1),xieta(2),'bo','LineWidth',2)

xk = isoparametricTransformationHighOrder(xieta,Xe,degree,nodesCoord)
