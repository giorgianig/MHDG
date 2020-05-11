% Postprocess routine

clear
% close all

global testcase Mesh

D = 0.1;
lscale = 1.901*1e-3;
% lscale = 1;
kT = 1;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
np = 1000;
r = 0.7501;

%% Solution with drift

%% Load mesh
meshName = 'West_h0.04_refCorn_P8_fixed.mat';
path2Mesh = '../../Meshes/Meshes_2D/WEST/';

%% Load solution
path2Sol = './Saves/';
solName = 'u_West_h04_refCorn_diff0.02_P8.mat';

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
Cline = Cline';

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
plotSolution(X,T,upar,refEl,20)
figure

% Circular line
Cline(:,1) = Cline(:,1)+shift;
Cline = Cline/lscale;
uplot = evalDGapproximationAtPoints(Cline,upar,X,T,refEl);
plot(Cline_theta,uplot,'k-'),grid on
ylabel('Mach')
xlabel('Angle')
xlim([min(Cline_theta), max(Cline_theta)]);
set(gca,'xtick',[pi/2 pi 3/2*pi 2*pi 5/2*pi] , ...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','-\pi'},'fontname','times new roman',...
    'fontsize',20)
set(gca,'ytick',-1.5:0.5:1.5,'fontname','times new roman',...
    'fontsize',20)
% readyforprintnew([8 6],24,[],[],1,[],[],'upar_on_separatrix_Drift')




