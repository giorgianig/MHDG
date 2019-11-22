% postprocess West Hole-NoHole
clear 
close all
global Mesh

lscale = 1.901*1e-3; 
npoints = 100;
Drift =1;
eV =200;
if eV == 50
    kT = 25;
    dfv = 2e-4;
elseif eV == 200
    kT = 100; 
    dfv = 8e-4;
end

Mesh.lscale = lscale;
%% Start
% if Drift
%     DriftStr = '_Drift';
% else
%     DriftStr = '';
% end



%% Find the separatrix
load ../../Meshes/Meshes_2D/WEST/WEST_wall.mat
load ../../Meshes/Meshes_2D/WEST/WEST_far_465.mat

%% Define another line
[Cline1,h] = contour(r2D,z2D,flux2D,[-0.86 -0.86]);
close
aux = sqrt((Cline1(1,:)-2.5).^2+(Cline1(2,:)+0.05).^2)< 0.75;
Cline1 = Cline1(:,aux);
df1 = diff(Cline1(1,:));
df2 = diff(Cline1(2,:));
dist = sqrt(df1.^2+df2.^2);
CurvCoo = [0,cumsum(dist)];



%% No hole solution
solpath = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/NoHole/';
if Drift
      solname = ['Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P8_Diff.38000E-01_YesDrift_', num2str(eV), 'eV.h5'];
else
    solname = ['Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P8_Diff.38000E-01_NoDrift_', num2str(eV), 'eV.h5'];
end
Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];
HDF5load([solpath,meshname])
HDF5load([solpath,solname])
u = u';
refEl = createReferenceElement(1,size(T,2));
Np = size(refEl.N,2);
Nhole.u = u; 
Nhole.X = X/lscale;
Nhole.T = T;

%% Hole solution
solpath = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/YesHole/';
if Drift
      solname = ['Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P8_Diff.38000E-01_YesDrift_', num2str(eV), 'eV.h5'];
else
    solname = ['Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P8_Diff.38000E-01_NoDrift_', num2str(eV), 'eV.h5'];
end
Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];
HDF5load([solpath,meshname])
HDF5load([solpath,solname])
u = u';
Yhole.u = u; 
Yhole.X = X/lscale;
Yhole.T = T;


%% Mapping
tol = 1e9;
x = Nhole.X(:,1); y = Nhole.X(:,2);
Nhole_xb = floor(1/3*sum(x(Nhole.T(:,1:3)),2)*tol);
Nhole_yb = floor(1/3*sum(y(Nhole.T(:,1:3)),2)*tol);
x = Yhole.X(:,1); y = Yhole.X(:,2);
Yhole_xb = floor(1/3*sum(x(Yhole.T(:,1:3)),2)*tol);
Yhole_yb = floor(1/3*sum(y(Yhole.T(:,1:3)),2)*tol);
[aux,map] = ismember([Yhole_xb,Yhole_yb],[Nhole_xb,Nhole_yb],'rows');

%% Work
uproj = zeros(size(Yhole.u));
for iel = 1:size(Yhole.T,1)
    
    ind_loc = (iel-1)*Np*2 + (1:2*Np);
    ind_map = (map(iel)-1)*Np*2 + (1:2*Np);
    uproj(ind_loc) = Nhole.u(ind_map);
end

% figure, plotSolution(Yhole.X,Yhole.T,Yhole.u(2:2:end)./Yhole.u(1:2:end-1)./sqrt(kT),refEl,5);caxis([-1.234762435076364   1.287260024619695])
figure, plotSolution(Yhole.X,Yhole.T,uproj(2:2:end)./uproj(1:2:end-1)./sqrt(kT)-...
                                      Yhole.u(2:2:end)./Yhole.u(1:2:end-1)./sqrt(kT) ,refEl,5);
% figure, plotSolution(Yhole.X,Yhole.T,uproj(2:2:end)./uproj(1:2:end-1)/dfv-...
%                                       Yhole.u(2:2:end)./Yhole.u(1:2:end-1)/dfv ,refEl,5);                                  
% figure, plotSolution(Yhole.X,Yhole.T,(uproj(1:2:end-1)-Yhole.u(1:2:end-1))./Yhole.u(1:2:end-1) ,refEl,5);                                  
% 
% figure, plotSolution(Yhole.X,Yhole.T,uproj(1:2:end-1),refEl,5);   
%                                caxis([-100 100])
caxis([-0.01 0.01])
% colorbar off
% figure, plotSolution(Nhole.X,Nhole.T,Nhole.u(1:2:end-1),refEl,5);   colorbar off

% figure, plotSolution(Yhole.X,Yhole.T,Yhole.u(1:2:end-1),refEl,5);   colorbar off

axis off      
% caxis([0 0.6])

readyforprintjpeg([8 6],24,[],[],[],[],[],'/home/giorgio/Dropbox/Conferences/PET_2017/Presentation/','DiffMachHoleNoHole_YesDrift_NEW')

% 
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,abs(Yhole.u(1:2:end-1)-uproj(1:2:end-1)),Yhole.X,Yhole.T,refEl);
% uplot(uplot==0) = NaN;
% figure, semilogy(line(:,1),uplot,'k'),grid on
% hold on
% uplot = evalDGapproximationAtPoints(line,Yhole.u(1:2:end-1),Yhole.X,Yhole.T,refEl);
% uplot(uplot==0) = NaN;
% semilogy(line(:,1),uplot,'k--'),grid on
% xlim([900 1600])
% ylim([1e-12 1])
% xlabel('Radial distance','fontname','times new roman','fontsize',20)
% ylabel('$\rho$','interpreter','latex','fontname','times new roman','fontsize',20)
% legend('Difference with vs without hole','With hole','Location','southeast')
%  readyforprintnew([8 6],24,[],[],[],[],[],'RhoProfile_HoleVsNohole_YesDrift')
% 
% 
% 
%  
%  
%  
 %% Extract profile on separatrix
%  fig=figure;
%  set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
% line = Cline1'/lscale;
% % uplot = evalDGapproximationAtPoints(line,Yhole.u(2:2:end)./Yhole.u(1:2:end-1)/sqrt(kT)-...
% %     uproj(2:2:end)./uproj(1:2:end-1)/sqrt(kT),Yhole.X,Yhole.T,refEl);
% uplot = evalDGapproximationAtPoints(line,Yhole.u(2:2:end)./Yhole.u(1:2:end-1)/dfv-...
%     uproj(2:2:end)./uproj(1:2:end-1)/dfv,Yhole.X,Yhole.T,refEl);
% uplot(uplot==0) = NaN;
% CurvCoo = CurvCoo(~isnan(uplot));
% CurvCoo = CurvCoo-min(CurvCoo);
% uplot = uplot(~isnan(uplot));
% yyaxis left, plot(CurvCoo/lscale,uplot,'k'),grid on
% ylabel('Difference with vs without hole','interpreter','latex','fontname','times new roman','fontsize',20)
% line = Cline1'/lscale;
% % uplot = evalDGapproximationAtPoints(line,Yhole.u(2:2:end)./Yhole.u(1:2:end-1)/sqrt(kT),Yhole.X,Yhole.T,refEl);
% uplot = evalDGapproximationAtPoints(line,Yhole.u(2:2:end)./Yhole.u(1:2:end-1)/dfv,Yhole.X,Yhole.T,refEl);
% uplot(uplot==0) = NaN;
% uplot = uplot(~isnan(uplot));
% yyaxis right, plot(CurvCoo/lscale,uplot,'k--'),grid on
% ylabel('$u_{\parallel}/ \|u_{\perp}\|$','interpreter','latex','fontname','times new roman','fontsize',20)
% legend('Difference with vs without hole','With hole','Location','south')
% xlabel('Distance on the curve','fontname','times new roman','fontsize',20)
% xlim([0 1800])
% 
% readyforprintnew([8 6],24,[],[],1,[],[],'UparProfile_HoleVsNohole_NoDrift_uperp')
% 
% 
