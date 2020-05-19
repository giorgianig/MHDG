% Postprocess circular case with limiter

close all
global Mesh

%% Physical parameters
lscale = 1.901*1e-3;
kT = 1;
testcase.n = 50;
axisym = 1;
fref = 7.28e6; % reference frequency
D = 0.02;
nref = 5;

%% Load mesh
path2Mesh = '../../Meshes/Meshes_2D/';


%%*************************************************************************
%%
%% Case with large limiter
%%
%%*************************************************************************

% meshName = 'Circel_LIM_small_RefCircle_h0.01_h0.1_P8.mat';
% path2Sol = './Saves/';
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
% 

% solName = 'u_Circel_LIM_small_RefCircle_h0.01_h0.1_P8.mat_Diffusion_0.01_tau1_tc60.mat';
% load([path2Sol solName])

% figure, plotMesh(X,T),readyforprintjpeg([8 6],16,[],[],[],[],[],'MeshCircSmallLim')
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'SL_dens_Diff1e-2_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off,caxis([-1 1])
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'SL_upar_Diff1e-2_NDr')
% 

% np = 100;
% line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% figure(100)
% semilogy(line(:,1),uplot,'r-'),grid on
% 
% lambda_n_fin=zeros(5,1);
% for i =1:6
%     Diff = 0.1/(10^(1/5))^(i-1);
%     solName = ['u_Circel_LIM_small_RefCircle_h0.01_h0.1_P8.mat_Diffusion_' num2str(Diff) '_tau1_tc60.mat'];
%     load([path2Sol solName])
%     line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
%     uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
%     ist = find(abs(line(:,1)-2.85/lscale) == min(abs(line(:,1)-2.85/lscale)));
%     ist = ist(1);
%     ien = find(abs(line(:,1)-2.86/lscale) == min(abs(line(:,1)-2.86/lscale)));
%     ien = ien(1);
%     lambda_n_fin(i) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));  
% end
% 
% 

%%*************************************************************************
%%
%% Case with infinitely small limiter
%%
%%*************************************************************************



% meshName = 'Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat';
% path2Sol = './Saves/';
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

% solName = 'u_Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_0.01_tau1_tc60.mat';
% load([path2Sol solName])
% 
% figure, plotMesh(X,T),readyforprintjpeg([8 6],16,[],[],[],[],[],'MeshCircInfThinLim')
% 
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_dens_Diff1e-2_NDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off,caxis([-1,1])
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82])
% readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_upar_Diff1e-2_NDr')
% 
% line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
% uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
% figure(100)
% hold on
% semilogy(line(:,1),uplot,'b-'),grid on
% ylabel('Density')
% xlabel('r')
% legend('Finite limiter','Infinitely thin limiter limiter','Location','Southwest')
% readyforprintnew([8 6],24,[],[],1,[],[],'Lambda_n_SmallLimVsInfSmallLimi')
% 
% lambda_n_ifin=zeros(5,1);
% for i =1:6
%     Diff = 0.1/(10^(1/5))^(i-1);
%     solName = ['u_Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_' num2str(Diff) '_tau1_tc60.mat'];
%     load([path2Sol solName])
%     line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
%     uplot = evalDGapproximationAtPoints(line,u0(1:2:end-1),X,T,refEl);
%     ist = find(abs(line(:,1)-2.85/lscale) == min(abs(line(:,1)-2.85/lscale)));
%     ist = ist(1);
%     ien = find(abs(line(:,1)-2.86/lscale) == min(abs(line(:,1)-2.86/lscale)));
%     ien = ien(1);
%     lambda_n_ifin(i) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));  
% end
% 
% figure(10),loglog(0.1./(10^(1/5)).^((1:6)-1),lambda_n_fin/lscale,'r--o',...
%                   0.1./(10^(1/5)).^((1:6)-1),lambda_n_ifin/lscale,'b--o');
% legend('Finite limiter','Infinitely thin limiter','location','northwest')
% ylabel('\lambda_n')
% xlabel('Diffusion')
% grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'LambdanVsDiff_SmallLimVsInfSmallLimi')






%% Case with drift
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
% path2Sol = './Saves/';
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
% 
% solName = 'u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat_Diffusion_0.0004_tau1_Drift_SC2_LocDiffPts2_NOCONV478.mat';
% load([path2Sol solName])
% % 
% %
% figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82],'fontname','times new roman')
% readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_smallSOL_dens_Diff4e-4_YDr')
% figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off
% colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82],'fontname','times new roman')
% readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_smallSOL_upar_Diff4e-4_YDr')
% 
% 
% close all


%% Case without drift
meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat';
path2Sol = './Saves/';
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

solName = 'u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat_Diffusion_0.0001_tau1_LocDiffPts2.mat';
load([path2Sol solName])
% 
%
figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref); axis off
colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82],'fontname','times new roman')
readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_smallSOL_dens_Diff4e-4_NDr')
figure,plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,nref); axis off,caxis([-1 1])
colorbar off, colorbar('location','east','position',[0.87,0.1,0.04,0.82],'fontname','times new roman')
readyforprintjpeg([8 6],16,[],[],[],[],[],'IL_smallSOL_upar_Diff4e-4_NDr')


close all