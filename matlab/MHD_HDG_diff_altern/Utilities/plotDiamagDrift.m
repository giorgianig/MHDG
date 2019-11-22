function plotDiamagDrift(X,T,varargin)

global Magnetic testcase Mesh

X = X*Mesh.lscale;
% if isempty(varargin)
%     fig = 1;
% else
%     fig = varargin{1};
% end

nElems = size(T,1);
Nv = size(T,2);

% initialize figure
% figure(fig),clf
% plotMesh(X,T)
% box on; axis equal tight;
% hold on

% loop in elements
for ielem = 1:nElems
    % index
    Te = T(ielem,:);
    if (testcase.n >= 50 && testcase.n<60)
        dx = Magnetic.dvxnodes(:,ielem);
        dy = Magnetic.dvynodes(:,ielem);
    else
        x = X(Te,1)/Mesh.lscale; y = X(Te,2)/Mesh.lscale;
       [~,~,dx,dy] = defineMagneticField([x,y]) ;
    end
    
    % plot
    quiver(X(Te,1),X(Te,2),dx,dy,100,'Color','b');
    hold on
end
axis equal; 
% axis([-0.1,1.1,-0.1,1.1]);
