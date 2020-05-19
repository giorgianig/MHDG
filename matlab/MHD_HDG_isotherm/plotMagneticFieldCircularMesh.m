% plot Magnetic field circular

clear
close 
path =  '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/';
mesh = 'Circel_LIM_small_RefCircle_h0.01_h0.1_P2.mat';
% mesh = 'Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P3.mat';
nlines = 20;


load([path mesh])

xmax = max(X(:,1));
xmin = min(X(:,1));
ymax = max(X(:,2));
ymin = min(X(:,2));
cent = 0.5*(xmax+xmin);

x = xmin*1.1:0.01:xmax*1.1;
y = ymin*1.1:0.01:ymax*1.1;

[XX,YY] = meshgrid(x,y);
mag = (sqrt((XX-cent).^2+YY.^2));

% surf(XX,YY,mag)

% aux = whos('Tb*');
% boundaryNames = {aux.name};
% clear aux
% nOfBound = numel(boundaryNames);
% for ibound = 1:nOfBound
%     name = boundaryNames{ibound};
%     eval(['aux_bound = Tb_' name(4:end) ';'])
%     for ifa = 1:size(aux_bound,1)
%         face = aux_bound(ifa,:);
%         plot(X(face,1),X(face,2),'k-','linewidth',2)
%         hold on
%     end
% end

plotMesh(X,T),hold on
axis equal,axis off

contour(XX,YY,mag,nlines,'k--','linewidth',0.5)


readyforprintnew([8 6],[],[],[],[],[],[],'Mesh_Circ_Magfield')


