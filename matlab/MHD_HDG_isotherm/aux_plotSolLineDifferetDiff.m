clear
% close all



lscale = 1.901*1e-3;
kT = 25;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
points = 1:7;
D = 0.1;
Diffvec = D*1./sqrt(10).^(points-1);
np = 100;
r = 0.3;
linestyle = {'k-', 'k--', 'k:','k-', 'k--', 'k:','k-', 'k--', 'k:'};
markstyle = {'ko-','k*--','ks:','k<-','k>--','kd:','kp-'};
hh = zeros(numel(points),1);
dth = 8*pi/180;
mks = 8; 
% theta0 = pi/2;
theta0 = -pi/2;
theta = theta0+linspace(0,2*pi,np)';
theta = theta(2:end-1);
%% ************************************************************************
%
%%                          Cases with Drift
%
%% ************************************************************************

disp('  ')
disp('Cases with drift')
disp('  ')

%% Diffusion 0.1 - 0.01

% Load mesh
meshName = 'Circle_ONION_1_P6.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

% Initialize
leg = {};
lambda_n_drift = zeros(size(points,2),1);
lambda_n_est_drift = lambda_n_drift;
leg_drift = cell(numel(points),1);

figure


disp('Diffusion: 0.1')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';
solName = 'u_Circle_ONION_1_P6_Diff.10000E+00.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{1}),grid on
hold on

ind = find(abs(theta-dth*1) == min(abs(theta-dth*1)));
hh(1)=plot(theta(ind),uplot(ind),markstyle{1},'markersize',mks);
leg{1} = num2str(Diffvec(1),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])











disp('Diffusion: 0.031623')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';
solName = 'u_Circle_ONION_1_P6_Diff.31623E-01.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{2}),grid on
hold on

ind = find(abs(theta-dth*2) == min(abs(theta-dth*2)));
hh(2)=plot(theta(ind),uplot(ind),markstyle{2},'markersize',mks);
leg{2} = num2str(Diffvec(2),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])













disp('Diffusion: 0.01')

meshName = 'Circle_ONION_2_P6.mat';
%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';
solName = 'u_Circle_ONION_2_P6_Diff.10000E-01.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{3}),grid on
hold on

ind = find(abs(theta-dth*3) == min(abs(theta-dth*3)));
hh(3)=plot(theta(ind),uplot(ind),markstyle{3},'markersize',mks);
leg{3} = num2str(Diffvec(3),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])






disp('Diffusion: 0.0031623')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-3/';
solName = 'u_Circle_ONION_2_P6_Diff.31623E-02_0076.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{4}),grid on
hold on

ind = find(abs(theta-dth*4) == min(abs(theta-dth*4)));
hh(4)=plot(theta(ind),uplot(ind),markstyle{4},'markersize',mks);
leg{4} = num2str(Diffvec(4),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])








disp('Diffusion: 0.001')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-3/';
solName = 'u_Circle_ONION_2_P6_Diff.10000E-02_0424.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{5}),grid on
hold on

ind = find(abs(theta-dth*5) == min(abs(theta-dth*5)));
hh(5)=plot(theta(ind),uplot(ind),markstyle{5},'markersize',mks);
leg{5} = num2str(Diffvec(5),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])











disp('Diffusion: 0.00031623')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-4/';
solName = 'u_Circle_ONION_2_P6_Diff.31623E-03_0448.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{6}),grid on
hold on

ind = find(abs(theta-dth*6) == min(abs(theta-dth*6)));
hh(6)=plot(theta(ind),uplot(ind),markstyle{6},'markersize',mks);
leg{6} = num2str(Diffvec(6),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])




meshName = 'Circle_ONION_2_P10.mat';
disp('Diffusion: 0.0001')

%% Load solution
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-4/';
solName = 'u_Circle_ONION_2_P10_Diff.10000E-03_0310.mat';
load([path2Mesh meshName])
load([path2Sol solName])

% Density and velocity
rho = u(1:2:end-1);
gamma = u(2:2:end);
upar = gamma./rho;

% Apply translation for axisymmetric case
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

[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar/sqrt(kT),X,T,refEl);
plot(theta,uplot,linestyle{7}),grid on
hold on

ind = find(abs(theta-dth*7) == min(abs(theta-dth*7)));
hh(7)=plot(theta(ind),uplot(ind),markstyle{7},'markersize',mks);
leg{7} = num2str(Diffvec(7),'%1.0e\n');
%     mTextBox = text(theta(ind),uplot(ind),num2str(Diffvec(idiff)));
%     set(mTextBox,'fontname','times new roman','BackgroundColor',[1 1 1])













ylabel('$M_{\|}$','interpreter','latex')
xlabel('\theta')
% xlim([min(theta), max(theta)]);
% set(gca,'xtick',[pi/2 pi 3/2*pi 2*pi 5/2*pi] , ...
%     'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'fontname','times new roman',...
%     'fontsize',20)

xlim([-pi/2 3/2*pi]);
set(gca,'xtick',[0 pi/2 pi 3/2*pi 2*pi]-pi/2, ...
   'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'},'fontname','times new roman',...
    'fontsize',20)


% set(gca,'ytick',-10:2:10,'fontname','times new roman',...
%     'fontsize',20)
legend(hh,leg,'location','best')
readyforprintnew([8 6],24,[],[],1,[],[],'upar_profile_r03_rev')
