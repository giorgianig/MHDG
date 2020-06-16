% Postprocess routine

clear
close all

global testcase Mesh

D = 0.1;
lscale = 1.901*1e-3;
% lscale = 1;
kT = 1;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
np = [100, 2000];
r = [0.75 0.7699];
linestyle = {'k-', 'k--'};
%% Solution with drift

%% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

%% Load solution
path2Sol = './Saves/';
solName = 'u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat_Diffusion_0.0001_tau1_Drift_SC2_LocDiffPts2.mat';



%% Working...
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u0(1:2:end-1);
gamma = u0(2:2:end);
upar = gamma./rho;


%% Apply translation if axisymmetric case
shift = 0;
if axisym && min(X(:,1))<0
    shift = - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    % apply translation in x to get min(X(:,1)) = 1
    X(:,1) = X(:,1) +shift;
end

X = X/lscale;

Mesh.X = X;
Mesh.T = T;
Mesh.lscale = lscale;


refEl = createReferenceElement(1,size(T,2),[]);

% figure
% plotSolution(X,T,upar,refEl,20)
% stop



figure
% Circular line

for iline = 1:2
    theta = pi/2+linspace(0,2*pi,np(iline))';
    [line_x, line_y] = pol2cart(theta,r(iline)*ones(size(theta)));    
    % line = ([line_x+shift, line_y] );
    line = ([line_x+shift, line_y] )/lscale;
    uplot = evalDGapproximationAtPoints(line,upar,X,T,refEl);
    plot(theta,uplot,linestyle{iline}),grid on
    hold on
end
ylabel('Mach')
xlabel('Angle')
xlim([min(theta), max(theta)]);
set(gca,'xtick',[pi/2 pi 3/2*pi 2*pi 5/2*pi] , ...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','-\pi'},'fontname','times new roman',...
    'fontsize',20)
set(gca,'ytick',-1.5:0.5:1.5,'fontname','times new roman',...
    'fontsize',20)
legend('On the separatrix','In the SOL','location','southeast')
% readyforprintnew([8 6],24,[],[],1,[],[],'upar_on_separatrix_Drift_NEW')







%% Solution without drift

%% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

%% Load solution
path2Sol = './Saves/';
solName = 'u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat_Diffusion_0.0001_tau1_LocDiffPts2.mat';

%% Working...
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u0(1:2:end-1);
gamma = u0(2:2:end);
upar = gamma./rho;

%% Apply translation if axisymmetric case
shift = 0;
if axisym && min(X(:,1))<0
    shift = - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    % apply translation in x to get min(X(:,1)) = 1
    X(:,1) = X(:,1) +shift;
end

X = X/lscale;

Mesh.X = X;
Mesh.T = T;
Mesh.lscale = lscale;

refEl = createReferenceElement(1,size(T,2),[]);
figure
% Circular line
for iline = 1:2
    theta = pi/2+linspace(0,2*pi,np(iline))';
    [line_x, line_y] = pol2cart(theta,r(iline)*ones(size(theta)));    
    % line = ([line_x+shift, line_y] );
    line = ([line_x+shift, line_y] )/lscale;
    uplot = evalDGapproximationAtPoints(line,upar,X,T,refEl);
    plot(theta,uplot,linestyle{iline}),grid on
    hold on
end
ylabel('Mach')
xlabel('Angle')
xlim([min(theta), max(theta)]);
set(gca,'xtick',[pi/2 pi 3/2*pi 2*pi 5/2*pi] , ...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','-\pi'},'fontname','times new roman',...
    'fontsize',20)
set(gca,'ytick',-1:0.2:1,'fontname','times new roman',...
    'fontsize',20)
legend('On the separatrix','In the SOL','location','northeast')

% readyforprintnew([8 6],24,[],[],1,[],[],'upar_on_separatrix_NoDrift_NEW')

