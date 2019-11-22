NZones=int32('/NZones'));

cmap=colormap('jet');
z=linspace(0,1,64);
z2=linspace(0,1,NZones);
% cmap2=[interp1(z,cmap(:,1),z2)',interp1(z,cmap(:,2),z2)',interp1(z,cmap(:,3),z2)'];

close all

for k=1:NZones
    Rcorner=hdf5read('mesh.h5',['/zone',num2str(k),'/Rcorner']);
    Zcorner=hdf5read('mesh.h5',['/zone',num2str(k),'/Zcorner']);
    Rgeom=hdf5read('mesh.h5',['/zone',num2str(k),'/Rgeom']);
    Zgeom=hdf5read('mesh.h5',['/zone',num2str(k),'/Zgeom']);
    figure(1)
    hold on
%     Rcorner = Rcorner(1:4:end,1:4:end);
%     Zcorner = Zcorner(1:4:end,1:4:end);
    
    plot(Rcorner,Zcorner,'-','Color','r');
    plot(Rcorner',Zcorner','-','Color','r');
%     plot(Rgeom(2:end-1,2:end-1),Zgeom(2:end-1,2:end-1),'.','Color','r');
end

% for k=1:NZones
%     figure(1)
%     hold on
%     Rgeom=hdf5read('mesh.h5',['/zone',num2str(k),'/Rgeom']);
%     Zgeom=hdf5read('mesh.h5',['/zone',num2str(k),'/Zgeom']);
%     [Nx,Nz]=size(Rgeom);
% %     text(Rgeom(floor(Nx/2),floor(Nz/2)),Zgeom(floor(Nx/2),floor(Nz/2)),num2str(k))
% end
% colormap(cmap2)
% colorbar
axis off
axis equal tight