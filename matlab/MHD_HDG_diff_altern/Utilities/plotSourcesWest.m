% plot geometry
home
clear
close
load ../../../Meshes/Meshes_2D/WEST/WEST_wall.mat
load ../../../Meshes/Meshes_2D/WEST/WEST_far_465.mat


%% No hole
% load ../../../Meshes/Meshes_2D/WEST/West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P2.mat
% flux2D(z2D<-0.7) = 1;

% Large source
% figure
% [C,h] = contourf(r2D,z2D,flux2D,[-0.90 -0.92]);
% hold on
% plot(Rwall,Zwall,'k-'); axis equal
% plotMesh(X,T);axis off
% pause(1)
% h.FacePrims(1).ColorData = uint8([255 0 0 255]'); 
% h.FacePrims(2).ColorData = uint8([255 255 255 255]');
% readyforprintnew([8 6],[],[],[],0.1,[],[],'figuresWest/NoHole_LS')

% Small source
% figure
% [C,h] = contourf(r2D,z2D,flux2D,[-1.03 -2] );
% hold on
% plot(Rwall,Zwall,'k-'); axis equal
% plotMesh(X,T);axis off
% pause(1)
% h.FacePrims(1).ColorData = uint8([255 0 0 255]'); 
% h.FacePrims(2).ColorData = uint8([255 255 255 255]');
% readyforprintnew([8 6],[],[],[],0.1,[],[],'figuresWest/NoHole_SS')

%% Small hole
% load ../../../Meshes/Meshes_2D/WEST/West_SmallHole_h0.04_refCorn_P6.mat
% 
% % Large source
% figure
% flux2D(z2D<-0.7) = 1;
% [C,h] = contourf(r2D,z2D,flux2D,[-0.90 -0.92]);
% hold on
% plot(Rwall,Zwall,'k-'); axis equal
% plotMesh(X,T);axis off
% pause(1)
% h.FacePrims(1).ColorData = uint8([255 0 0 255]'); 
% h.FacePrims(2).ColorData = uint8([255 255 255 255]');
% readyforprintnew([8 6],[],[],[],0.1,[],[],'figuresWest/SmallHole_LS')

% Small source
% figure
% flux2D(z2D<-0.7) = 1;
% [C,h] = contourf(r2D,z2D,flux2D,[-0.98 -1]);
% hold on
% plot(Rwall,Zwall,'k-'); axis equal
% plotMesh(X,T);axis off
% pause(1)
% h.FacePrims(1).ColorData = uint8([255 0 0 255]'); 
% h.FacePrims(2).ColorData = uint8([255 255 255 255]');
% readyforprintnew([8 6],[],[],[],0.1,[],[],'figuresWest/SmallHole_SS')


% Large hole
load ../../../Meshes/Meshes_2D/WEST/West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P2.mat

figure
flux2D(z2D<-0.7) = 1;
flux2D(z2D>0.64) = 1;

[C,h] = contourf(r2D,z2D,flux2D,[-0.88 -0.90]);
hold on
plot(Rwall,Zwall,'k-'); axis equal
plotMesh(X,T);axis off
pause(1)
h.FacePrims(1).ColorData = uint8([255 0 0 255]'); 
h.FacePrims(2).ColorData = uint8([255 255 255 255]');
% readyforprintnew([8 6],[],[],[],0.1,[],[],'figuresWest/LargeHole_LS')
readyforprintjpeg([8 6],[],[],[],0.1,[],[],'MeshWestNoHoleSource')
