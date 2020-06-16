% Postprocess routine

clear
close all

global testcase Mesh

%% Physical parameters
lscale = 1.901*1e-3;
kT = 1;
testcase.n = 50;
axisym = 1;
fref = 7.28e6; % reference frequency
D = 0.02;
nref = 5;
npoints = 100;
%% Find the separatrix
load ../../Meshes/Meshes_2D/WEST/WEST_wall.mat
load ../../Meshes/Meshes_2D/WEST/WEST_far_465.mat
[Cline,h] = contour(r2D,z2D,flux2D,[-0.8662 -0.8662]);
close
aux = sqrt((Cline(1,:)-2.5).^2+Cline(2,:).^2)<sqrt((2.5-2.2722)^2+0.6082^2);
Cline = Cline(:,aux);
Centx = mean(Rwall);
Centy = mean(Zwall);
Cline_theta = pi/2+cart2pol(Cline(1,:)-Centx,Cline(2,:)-Centy);
Cline_theta(Cline_theta>pi)= Cline_theta(Cline_theta>pi)-2*pi;
[Cline_theta,ind] = sort(Cline_theta,'ascend');
Cline = Cline(:,ind);

%% Load mesh
path2Mesh = '../../Meshes/Meshes_2D/WEST/';
path2Sol = './Saves/';

%%*************************************************************************
%%
%% Large hole
%%
%%*************************************************************************
% meshName = 'West_h0.04_refCorn_P8_fixed.mat';

% load([path2Mesh meshName])
% 
% % Apply translation if axisymmetric case
% if axisym && min(X(:,1))<0
%     % apply translation in x to get min(X(:,1)) = 1
%    X(:,1) = X(:,1) - min(X(:,1))+1;
% end
% refEl = createReferenceElement(1,size(T,2),[]);
% X = X/lscale;
% 
% Mesh.X = X;
% Mesh.T = T;
% Mesh.lscale = lscale;

% CASE WITH Dirichlet BC
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc50.mat';
% load([path2Sol solName])
% % 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off,colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_LH_NS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_NS_NDr_NoCb_AxisAuto')


% Colorbar
% figure, colorbar('location','north','position',[0.1 0.1 0.85 0.04],'fontname','times new roman','fontsize',16)
% colormap('jet'), caxis([-1.289 1.317]),axis off,readyforprintjpeg([8 6],20,[],[],[],[],[],'Cbar_Upar')

% stop

% CASE WITHOUT DRIFT
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_LH_LS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_NDr_NoCb_AxisAuto_Super')


% Colorbar
% figure, colorbar('location','north','position',[0.1 0.1 0.85 0.04],'fontname','times new roman','fontsize',16)
% colormap('jet'), caxis([-1.289 1.317]),axis off,readyforprintjpeg([8 6],20,[],[],[],[],[],'Cbar_Upar')

% Plot with separatrix
% plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off, colorbar off
% hold on
% plot(Cline(1,:),Cline(2,:),'k--')
% readyforprintjpeg([8 6],16,[],[],1,[],[],'Upar_LH_LS_NDr_NoCb_AxisAuto_Super_Separatrix')

% Plot on the large hole mesh the difference of the solution between large hole
% and no hole
% solNameLargeHole = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51.mat';
% solNameNoHole    = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_PROJECTEDONLARGEHOLE.mat';
% solLH = load([path2Sol solNameLargeHole]);
% solNH = load([path2Sol solNameNoHole]);
% figure,plotSolution(X,T,(solLH.u0(1:2:end-1)-solNH.upr(1:2:end-1))./solNH.upr(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.01])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_LH_NH_LS_NDr_RelativeNH_scale001')
% figure,plotSolution(X,T,solLH.u0(1:2:end-1)-solNH.upr(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.01])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_LH_NH_LS_NDr_scale001')
% figure,plotSolution(X,T,(solLH.u0(2:2:end)./solLH.u0(1:2:end-1)-solNH.upr(2:2:end)./solNH.upr(1:2:end-1))... 
%                                        ./(solNH.upr(2:2:end)./solNH.upr(1:2:end-1)),refEl,nref); axis off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'UparDiff_LH_NH_LS_NDr_RelativeNH')
% figure,plotSolution(X,T,(solLH.u0(2:2:end)./solLH.u0(1:2:end-1)-solNH.upr(2:2:end)./solNH.upr(1:2:end-1))... 
%                                        ,refEl,nref); axis off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'UparDiff_LH_NH_LS_NDr')


% CASE WITH DRIFT
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_LH_LS_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr_NoCb_AxisAuto_Super')

% Plot on the large hole mesh the difference of the solution between large hole
% and no hole
% solNameLargeHole = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% solNameNoHole    = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift_PROJECTEDONLARGEHOLE.mat';
% solLH = load([path2Sol solNameLargeHole]);
% solNH = load([path2Sol solNameNoHole]);
% figure,plotSolution(X,T,(solLH.u0(1:2:end-1)-solNH.upr(1:2:end-1))./solNH.upr(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.01])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_LH_NH_LS_YDr_RelativeNH_scale001')
% figure,plotSolution(X,T,(solLH.u0(1:2:end-1)-solNH.upr(1:2:end-1)),refEl,nref); axis off,caxis([-0.01 0.01])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_LH_NH_LS_YDr_scale001')
% figure,plotSolution(X,T,(solLH.u0(2:2:end)./solLH.u0(1:2:end-1)-solNH.upr(2:2:end)./solNH.upr(1:2:end-1))... 
%                                        ./(solNH.upr(2:2:end)./solNH.upr(1:2:end-1)),refEl,nref); axis off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'UparDiff_LH_NH_LS_YDr_RelativeNH')
% figure,plotSolution(X,T,(solLH.u0(2:2:end)./solLH.u0(1:2:end-1)-solNH.upr(2:2:end)./solNH.upr(1:2:end-1))... 
%                                        ,refEl,nref); axis off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'UparDiff_LH_NH_LS_YDr')

% stop
% DIFFERECES WITH/WITHOUT DRIFT
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51.mat';
% und = load([path2Sol solName]);
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% uyd = load([path2Sol solName]);
% figure,plotSolution(X,T,und.u0(1:2:end-1)-uyd.u0(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.02])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_LH_LS_YNDr_NCb')
% figure, colorbar('location','north','position',[0.1 0.1 0.85 0.04],'fontname','times new roman','fontsize',16)
% colormap('jet'), caxis([-0.01 0.02]),axis off,readyforprintjpeg([8 6],20,[],[],[],[],[],'Cbar_DensDiff')
% % figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% % readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')


% plotParallelVelocityVector(X,T,u(2:2:end)./u(1:2:end-1),10)

% ************************* Lambda n **************************************
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% 
% figure(100)
% legend_rhoprof = {};
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),uplot,'r-','Markersize',16),grid on
% uplotref = uplot;
% legend_rhoprof{1} = 'Large hole - no drift - REFERENCE';

% figure(200)
% legend_uprof = {};
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'r-','Markersize',16),grid on
% legend_uprof{1} = 'Large hole - no drift';

% % % % 
% % % % line(uplot==0) = [];
% % % % uplot(uplot==0) = [];
% % % % 
% % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % Diff = D*fref*(lscale)^2;
% % % % vref = lscale*fref;
% % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % leg{1} = 'Large Hole';

% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'r--','Markersize',16),grid on
% legend_rhoprof{2} = 'Large hole - drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'r--','Markersize',16),grid on
% legend_uprof{2} = 'Large hole - drift';

%%*************************************************************************
%%
%% Small hole
%%
%%*************************************************************************

% meshName = 'West_SmallHole_h0.04_refCorn_P8.mat';
% load([path2Mesh meshName])
% 
% % Apply translation if axisymmetric case
% if axisym && min(X(:,1))<0
%     % apply translation in x to get min(X(:,1)) = 1
%    X(:,1) = X(:,1) - min(X(:,1))+1;
% end
% 
% X = X/lscale;
% refEl = createReferenceElement(1,size(T,2),[]);
% 
% Mesh.X = X;
% Mesh.T = T;
% Mesh.lscale = lscale;

% CASE WITHOUT DRIFT  - Source on the large ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_SH_LS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_SH_LS_NDr_NoCb_AxisAuto_Super')


% CASE WITHOUT DRIFT  - Source on the small ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_SH_SS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_SH_SS_NDr_NoCb_AxisAuto_Super')


% CASE WITH DRIFT  - Source on the large ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_SH_LS_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off 
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_SH_LS_YDr_NoCb_AxisAuto_Super')


% CASE WITH DRIFT  - Source on the small ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54_Drift.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_SH_SS_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , colorbar off 
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_SH_SS_YDr_NoCb_AxisAuto_Super')

% DIFFERECES WITH/WITHOUT DRIFT - Source on the large ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% und = load([path2Sol solName]);
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% uyd = load([path2Sol solName]);
% figure,plotSolution(X,T,und.u0(1:2:end-1)-uyd.u0(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.02])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_SH_LS_YNDr_NCb')
% % figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% % readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')

% DIFFERECES WITH/WITHOUT DRIFT - Source on the small ring
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54.mat';
% und = load([path2Sol solName]);
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54_Drift.mat';
% uyd = load([path2Sol solName]);
% figure,plotSolution(X,T,und.u0(1:2:end-1)-uyd.u0(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.02])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_SH_SS_YNDr_NCb')
% % figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% % readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')


% plotParallelVelocityVector(X,T,u(2:2:end)./u(1:2:end-1),10)


% ************************* Lambda n **************************************
% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% figure(100)
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'b-','Markersize',16),grid on
% legend_rhoprof{3} = 'Small hole - large source - no drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'b-','Markersize',16),grid on
% legend_uprof{3} = 'Small hole - large source - no drift';

% 
% % % % line(uplot==0) = [];
% % % % uplot(uplot==0) = [];
% % % % 
% % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % Diff = D*fref*(lscale)^2;
% % % % vref = lscale*fref;
% % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % leg{1} = 'Large Hole';

% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'b--','Markersize',16),grid on
% legend_rhoprof{4} = 'Small hole - large source - drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'b--','Markersize',16),grid on
% legend_uprof{4} = 'Small hole - large source - drift';


% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'g-','Markersize',16),grid on
% legend_rhoprof{5} = 'Small hole - small source - no drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'g-','Markersize',16),grid on
% legend_uprof{5} = 'Small hole - small source - no drift';

% % % % % 
% % % % line(uplot==0) = [];
% % % % uplot(uplot==0) = [];
% % % % 
% % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % Diff = D*fref*(lscale)^2;
% % % % vref = lscale*fref;
% % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % leg{1} = 'Large Hole';

% solName = 'u_West_SmallHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc54_Drift.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'g--','Markersize',16),grid on
% legend_rhoprof{6} = 'Small hole - small source - drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'g--','Markersize',16),grid on
% legend_uprof{6} = 'Small hole - small source - drift';


% % % % % % line(uplot==0) = [];
% % % % % % uplot(uplot==0) = [];
% % % % % % 
% % % % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % % % Diff = D*fref*(lscale)^2;
% % % % % % vref = lscale*fref;
% % % % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % % % leg{2} = 'Small Hole';

%%*************************************************************************
%%
%% No hole
%%
%%*************************************************************************

meshName = 'West_NoHole_h0.04_refCorn_P8.mat';
load([path2Mesh meshName])

% Apply translation if axisymmetric case
if axisym && min(X(:,1))<0
    % apply translation in x to get min(X(:,1)) = 1
   X(:,1) = X(:,1) - min(X(:,1))+1;
end
refEl = createReferenceElement(1,size(T,2),[]);

X = X/lscale;

Mesh.X = X;
Mesh.T = T;
Mesh.lscale = lscale;

% CASE WITHOUT DRIFT  - Source on the large ring
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% % 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_NH_LS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_NH_LS_NDr_NoCb_AxisAuto_Super')

% CASE WITHOUT DRIFT  - Source in the center
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_NH_SS_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_NH_SS_NDr_NoCb_AxisAuto_Super')


% CASE WITH DRIFT  - Source on the large ring
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_NH_LS_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_NH_LS_YDr_NoCb_AxisAuto_Super')

% CASE WITH DRIFT  - Source in the center
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55_Drift.mat';
% load([path2Sol solName])
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off, caxis([0,2])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Dens_NH_SS_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off, colorbar off
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_NH_SS_YDr_NoCb_AxisAuto_Super')


% DIFFERECES WITH/WITHOUT DRIFT - Source on the large ring
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% und = load([path2Sol solName]);
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% uyd = load([path2Sol solName]);
% figure,plotSolution(X,T,und.u0(1:2:end-1)-uyd.u0(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.02])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_NH_LS_YNDr_NCb')
% % figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% % readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')

% DIFFERECES WITH/WITHOUT DRIFT - Source on the center
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55.mat';
% und = load([path2Sol solName]);
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55_Drift.mat';
% uyd = load([path2Sol solName]);
% figure,plotSolution(X,T,und.u0(1:2:end-1)-uyd.u0(1:2:end-1),refEl,nref); axis off,caxis([-0.01 0.02])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_NH_SS_YNDr_NCb')
% % figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% % readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')

% DIFFERECES WITH DRIFT large source/source on the center
[Cline1,h] = contour(r2D,z2D,flux2D,[-0.90 -0.90]);
aux = sqrt((Cline1(1,:)-2.5).^2+Cline1(2,:).^2)<0.5;
Cline1 = Cline1(:,aux);
[Cline2,h] = contour(r2D,z2D,flux2D,[-0.92 -0.92]);
Cline2 = Cline2(:,2:end);
[Cline3,h] = contour(r2D,z2D,flux2D,[-1.03 -1.03]);
Cline3 = Cline3(:,2:end);
close

solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55_Drift.mat';
uss = load([path2Sol solName]);
solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
uls = load([path2Sol solName]);
rhoss= uss.u0(1:2:end-1); rhols = uls.u0(1:2:end-1);
upss= uss.u0(2:2:end)./uss.u0(1:2:end-1); upls = uls.u0(2:2:end)./uls.u0(1:2:end-1);

% figure,plotSolution(X,T,(rhoss-rhols)./rhoss,refEl,nref); axis off
% hold on, plot(Cline1(1,:),Cline1(2,:),'k--')
% hold on, plot(Cline2(1,:),Cline2(2,:),'k--')
% hold on, plot(Cline3(1,:),Cline3(2,:),'k--')
% stop
figure,plotSolution(X,T,(upss-upls)./upss,refEl,nref); axis off
hold on, plot(Cline1(1,:),Cline1(2,:),'k--')
hold on, plot(Cline2(1,:),Cline2(2,:),'k--')
hold on, plot(Cline3(1,:),Cline3(2,:),'k--')

% readyforprintjpeg([8 6],16,[],[],[],[],[],'DensDiff_NH_LS_YNDr_NCb')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref,0,0,1); axis off , caxis([-0.1,0.1])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'Upar_LH_LS_YDr')






% plotParallelVelocityVector(X,T,u(2:2:end)./u(1:2:end-1),10)

% % ************************* Lambda n **************************************

% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% figure(100)
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'k-','Markersize',16),grid on
% legend_rhoprof{7} = 'Small hole - large source - no drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'k-','Markersize',16),grid on
% legend_uprof{7} = 'Small hole - large source - no drift';

% % % % % % line(uplot==0) = [];
% % % % % % uplot(uplot==0) = [];
% % % % % % 
% % % % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % % % Diff = D*fref*(lscale)^2;
% % % % % % vref = lscale*fref;
% % % % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % % % leg{1} = 'Large Hole';


% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'k--','Markersize',16),grid on
% legend_rhoprof{8} = 'Small hole - large source - drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'k--','Markersize',16),grid on
% legend_uprof{8} = 'Small hole - large source - drift';
 

% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55.mat';
% load([path2Sol solName])
% figure(100)
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'m-','Markersize',16),grid on
% legend_rhoprof{9} = 'No hole - small source - no drift - Difference';

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'m-','Markersize',16),grid on
% legend_uprof{9} = 'No hole - small source - no drift';
% 
% % % % % % line(uplot==0) = [];
% % % % % % uplot(uplot==0) = [];
% % % % % % 
% % % % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % % % Diff = D*fref*(lscale)^2;
% % % % % % vref = lscale*fref;
% % % % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % % % leg{1} = 'Large Hole';

% 
% solName = 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc55_Drift.mat';
% load([path2Sol solName])
% figure(100)
% hold on
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% semilogy(line(:,1),abs(uplot-uplotref),'m--','Markersize',16),grid on
% legend_rhoprof{10} = 'No hole - small source - drift - Difference';
% legend(legend_rhoprof,'location','southeast')
% ylabel('Density')
% xlabel('r')
% readyforprintnew([12 9],24,[],[],1,[],[],'figuresWest/denslineWest')
% 
% pause(1)

% figure(200)
% hold on
% uplot = evalDGapproximationAtPoints(Cline'/lscale,u0(2:2:end)./u0(1:2:end-1),X,T,refEl);
% plot(Cline_theta,uplot,'m--','Markersize',16),grid on
% legend_uprof{10} = 'No hole - small source - drift';
% legend(legend_uprof,'location','northwest')
% set(gca,'xtick',[-pi -pi/2 0 pi/2 pi] , ...
%     'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','-\pi'},'fontname','times new roman',...
%     'fontsize',20)
% set(gca,'ytick',-1.5:0.1:1.5,'fontname','times new roman',...
%     'fontsize',20)
% ylabel('Mach')
% xlabel('Angle')
% xlim([min(Cline_theta)-0.01, max(Cline_theta)+0.01]);
% readyforprintnew([12 9],24,[],[],1,[],[],'figuresWest/uparsepar')

% % % % % % % line(uplot==0) = [];
% % % % % % % uplot(uplot==0) = [];
% % % % % % % 
% % % % % % % ist = find(abs(line(:,1)-2.8/lscale) == min(abs(line(:,1)-2.8/lscale)));
% % % % % % % ien = find(abs(line(:,1)-2.9/lscale) == min(abs(line(:,1)-2.9/lscale)));
% % % % % % % lambda_n = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
% % % % % % % Diff = D*fref*(lscale)^2;
% % % % % % % vref = lscale*fref;
% % % % % % % lambda_n_est = sqrt(3/4*43*Diff/vref);
% % % % % % % disp(['Lambda_n from graphic: ' num2str(lambda_n)])
% % % % % % % disp(['Lambda_n estimated: ' num2str(lambda_n_est)])
% % % % % % % leg{3} = 'No Hole';
% % % % % % % 
% % % % % % % legend(leg)
% % % % % % % xlabel('Radius')
% % % % % % % ylabel('Density')



% close all