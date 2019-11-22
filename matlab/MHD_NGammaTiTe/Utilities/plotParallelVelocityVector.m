function plotParallelVelocityVector(X,T,u,varargin)

global Magnetic testcase Mesh

X = X*Mesh.lscale;
if isempty(varargin)
    fig = 1;
else
    fig = varargin{1};
end

nElems = size(T,1);
Nv = size(T,2);

% initialize figure
figure(fig),clf
plotMesh(X,T)
box on; axis equal tight;
hold on

% loop in elements
for ielem = 1:nElems
    % index
    ind = (ielem-1)*Nv + (1:Nv);
    Te = T(ielem,:);
    if (testcase.n >= 50 && testcase.n<60)
        bx = Magnetic.bxnodes(:,ielem);
        by = Magnetic.bynodes(:,ielem);
    else
        x = X(Te,1)/Mesh.lscale; y = X(Te,2)/Mesh.lscale;
       b = defineMagneticField([x,y]) ;
       bx = b(:,1); by = b(:,2);
    end
    
    % plot
    quiver(X(Te,1),X(Te,2),u(ind).*bx,u(ind).*by,0,'Color','b');
    hold on
end
axis equal; 
% axis([-0.1,1.1,-0.1,1.1]);
