% plot Solution of moving equilibrium
clear all
close all
global Mesh


kT = 25;
neq = 2;
lscale = 1.901*1e-3; 
zoom = 0;
simu = 'slow';
% simu = 'fast';

mks = 20;
if strcmpi(simu,'fast')
    tref = 5e-5;
else
    tref = 5e-4;
end
if zoom
    nt = 90;
else
    if strcmpi(simu,'fast')
        nt = 2836;
    else
        nt =753;
    end
end
% solpath = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/NoHole/';
solpath = '/home/giorgio/Desktop/MHDG_Thomas/test/Evolv_Giorgio/Solutions_dte-1/';
meshpath = '/home/giorgio/Desktop/MHDG_Thomas/test/Evolv_Giorgio/';
mfieldpath =  '/home/giorgio/Desktop/MHDG_Thomas/test/Evolv_Giorgio/';
solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P3_Diff.38000E-01';
mfieldname = 'Evolving_equilibrium';
load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/Evolving_equilibrium_35132els_10ts_P3/Lim2Div/psi_int_10.mat
savepath  = '/home/giorgio/Dropbox/Conferences/PET_2017/figures/';
% start
Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];

HDF5load([meshpath,meshname])
refEl = createReferenceElement(1,size(T,2));
Nv = size(T,2);
Mesh.Tb = Tb;

% % Determine colorbar
% max_rho = 0;
% min_up = 100;
% max_up = -100;
% max_psi = -100;
% min_psi  = 100;
% for it =1:nt
%     solname_full = [solname,'_',num2str(it,'%04d'),'.h5'];
%     HDF5load([solpath,solname_full])
%     u = u';
%     max_rho = max(max_rho,max(u(1:2:end-1)));
%     max_up = max(max_up,max(u(2:2:end)./u(1:2:end-1)/sqrt(kT)));
%     min_up = min(min_up,min(u(2:2:end)./u(1:2:end-1)/sqrt(kT)));
%     max_psi = max(max_psi,max(psi(:,it)));
%     min_psi = min(min_psi,min(psi(:,it)));
% 
% end
% 
% disp(['Max rho:    ' num2str(max_rho)])
% disp(['Max Mach: ' num2str(max_up)])
% disp(['Min Mach:  ' num2str(min_up)])
% disp(['Max Psi:  ' num2str(max_psi)])
% disp(['Min Psi:  ' num2str(min_psi)])
% stop

[intf,extf] = GetFaces(T(:,1:3));

%% Find left face
% fcl = zeros(size(extf,1),1);
% ind = 0;
% indnod_l = zeros(size(extf,1),size(refEl.faceNodes,2));
% indnod_l_u = zeros(size(extf,1),size(refEl.faceNodes,2));
% for ifa = 1:size(extf,1)
%     nodes_loc = refEl.faceNodes(extf(ifa,2),:);
%     nodes = T(extf(ifa,1),refEl.faceNodes(extf(ifa,2),:));
%     if all(X(nodes,1)<1.865)
%         ind = ind+1;
%         fcl(ind) = ifa;
%         indnod_l(ind,:) = nodes;
%         indnod_l_u(ind,:) =  (extf(ifa,1)-1)*Nv + nodes_loc;
%     end
% end
% fcl = fcl(1:ind,1);
% indnod_l = col(transpose(indnod_l(1:ind,:)));
% indnod_l_u = col(transpose(indnod_l_u(1:ind,:)));
% [~,a] = sort(X(indnod_l,2));
% indnod_l = indnod_l(a);
% indnod_l_u = indnod_l_u(a,:);
% sl = sqrt(diff(X(indnod_l,1)).^2 + diff(X(indnod_l,2)).^2);
% sl = [0;cumsum(sl)];
% indnod_l_u = col(transpose(cat(2,indnod_l_u,2*indnod_l_u)));


for it =1:1
    solname_full = [solname,'_',num2str(it,'%04d'),'.h5'];
    HDF5load([solpath,solname_full])
    u = u';
%      
%     %% plot Psi
% %     h = figure('visible','off');  plotSolution(X,T,psi(:,it),refEl,5,1), axis off,caxis([ 1.0644  1.2286]),colorbar off
% %     readyforprintjpeg([8 6],24,[],[],[],h,[],savepath,['psi_cont_' num2str(it,'%04d')])
% %     close(h)
% % set(h,'visible','on')
% % 
    %% plot density
%      h = figure('visible','off'); plotSolution(X/lscale,T,u(1:2:end-1),refEl,5);axis off,caxis([0 2.1]),%colorbar off
%      readyforprintjpeg([8 6],24,[],[],[],h,[],savepath,['rho_slow_' num2str(it,'%04d')])
%      close(h)
% % %     
% %     %% plot Mach
%      h = figure('visible','off'); 
     plotSolution(X/lscale,T,u(2:2:end)./u(1:2:end-1)./sqrt(kT),refEl,5);axis off,caxis([-2   2]);%colorbar off
% %      readyforprintjpeg([8 6],24,[],[],[],h,[],savepath,['Mach_slow_' num2str(it,'%04d')])
% %      close(h)
% %      
%     %% plot shock-capturing coefficient
% %      h = figure('visible','off'); plotSolution(X/lscale,T,col(scdiff_nodes'),refEl,5);axis off,caxis([0   4e-2]);colorbar off
% %      readyforprintjpeg([8 6],24,[],[],[],h,[],solpath,['SC_coeff_' num2str(it,'%04d')])
% %      close(h)
% 
%   %% plot Outflux
% %   mfieldname_full = [mfieldname,'_',num2str(it,'%04d'),'.h5'];
% %   HDF5load([mfieldpath,mfieldname_full]) 
% %  h = figure('visible','off');
% %  plotSolutionOnBoundary(u(2:2:end),X,T,extf,refEl,[Br',Bz',Bt'],5,50,1)
% %  xlim([1   4])
% %  ylim([-1.5    1.5])
% %  readyforprintjpeg([8 6],24,[],[],[],h,[],solpath,['Outflux_' num2str(it,'%04d')])
% %  close(h)
% 
%  %% 2D plot of Gamma vs time
%    mfieldname_full = [mfieldname,'_',num2str(it,'%04d'),'.h5'];
%   HDF5load([mfieldpath,mfieldname_full]) 
%   [uint,xyint,indpt] = analyseBoundaryWest(u(2:2:end),X,T,extf,refEl,[Br',Bz',Bt'],5,1);
%   pl2d(it,:) = uint;
% 
end
% if strcmpi(simu,'fast')
%     load pl2d_fast.mat
% else
%     load pl2d_slow.mat
% end
% 
% 
% pl2d(pl2d<0) = 0;
% xdif = diff(xyint(:,1)); ydif = diff(xyint(:,2));
% xy = [0;cumsum( sqrt(xdif.^2+ydif.^2)   )];
% if strcmpi(simu,'fast')
%     if zoom
%        timevec = (1:nt)*tref;
%     else
%         timevec = [(1:2000)*tref, 2000*tref+(1:(nt-2000))*tref*10];
%     end
% else
%     timevec = (1:nt)*tref;
% end
% if timevec(end)>0.35
%     nt = find(abs(timevec-0.35) == min(abs(timevec-0.35)));
%     timevec = timevec(1:nt);
% end
% [xx,tt] = meshgrid(xy,timevec);
% contourf(xx,tt,pl2d(1:nt,:)),colormap('jet'),colorbar,caxis([0,5e-3])
% hold on
% plot(transpose(xy(indpt)*ones(1,nt)), (timevec'*ones(1,size(indpt,1))),'r--','linewidth',2)
% set(gca,'xtick',[0.25, 1.2,2])
% set(gca,'xticklabel',{'upper divertor','left limiter','lower divertor'},'fontname','times new roman')
% ylabel('time [s]','fontname','times new roman')
% 
% if zoom
%     if strcmpi(simu,'fast')
%         readyforprintjpeg([8 6],mks,[],[],[],[],[],savepath,'Outfluxvstime_fast_scale')
%     else
%         readyforprintjpeg([8 6],mks,[],[],[],[],[],savepath,'Outfluxvstime_slow_scale')
%     end    
% else
%     if strcmpi(simu,'fast')
%         readyforprintjpeg([8 6],mks,[],[],[],[],[],savepath,'Outfluxvstime_fast_scale_long')
%     else
%         readyforprintjpeg([8 6],mks,[],[],[],[],[],savepath,'Outfluxvstime_slow_scale_long')
%     end
% end
