% plot evolving equilibrium

global Mesh
Mesh.lscale = 1;
clear 
close all
path = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/Evolving_equilibrium_35132els_10ts_P3/Lim2Div/';
load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/WEST_far_465.mat
load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/WEST_wall.mat
load([path,'Coordinates.mat'])
load([path,'Elements.mat'])
load([path,'fpsi.mat'])
load([path,'psir.mat'])
load([path,'psiz.mat'])
load([path,'psi.mat'])


% compute axis
minx = min(Rwall);
maxx = max(Rwall);
miny = min(Zwall);
maxy = max(Zwall);

nt = size(fpsi,2);
refEl = createReferenceElement(1,size(T,2));

% compute bx, by 
r = X(:,1);
Bmod = sqrt(psir.^2+psiz.^2+fpsi.^2);
bx = -psiz./Bmod;
by =  psir./Bmod;

Btor = fpsi./repmat(r,[1,10]);
Bpol = sqrt( (psiz./repmat(r,[1,10])).^2 + (psir./repmat(r,[1,10])).^2);



% quiver(X(:,1),X(:,2),bx(:,1),by(:,1))

for it =[1,4,7,10]
%     ax = find(Bpol(:,it)==min(Bpol(:,it)));
%     Btor(ax,it)
    figure
    plot(Rwall,Zwall,'k');axis equal
    hold on
% %     figure, plotSolution(X,T,psi(:,it),refEl,5)
    [CC,hh]= plotSolution(X,T,psi(:,it),refEl,5,1,0,0) ;
% % quiver(X(:,1),X(:,2),bx(:,it),by(:,it)),axis equal
    axis([minx maxx miny maxy])
    axis off
end


% Comparison with static equilibrium
% Bmod_stat = sqrt(Br2D.^2+Bz2D.^2+Bphi2D.^2);
% bx_stat = Br2D./Bmod_stat;
% by_stat = Bz2D./Bmod_stat;
% figure, plot(Rwall,Zwall,'k');axis equal,hold on, contour(r2D,z2D,flux2D)