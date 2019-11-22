function plotVelocity(X,T,u,varargin)

neq = 3;
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
    indr = (ielem-1)*neq*Nv + (1:neq:neq*Nv-2);
    indu = (ielem-1)*neq*Nv + (2:neq:neq*Nv-1);
    indv = (ielem-1)*neq*Nv + (neq:neq:neq*Nv);
    Te = T(ielem,:);
    
    % plot
    quiver(X(Te,1),X(Te,2),u(indu)./u(indr),u(indv)./u(indr),0.5);
    hold on
end
axis equal; 
% axis([-0.1,1.1,-0.1,1.1]);