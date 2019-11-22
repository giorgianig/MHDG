clear
close all

global Mesh
shockcapt = 2;
neq=2;
lscale = 1.901*1e-3;
Mesh.lscale = lscale;

meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';
load([path2Mesh meshName])

% apply translation in x to get min(X(:,1)) = 1
shift = - 0.5*(min(X(:,1))+max(X(:,1)))+2;
% apply translation in x to get min(X(:,1)) = 1
X(:,1) = X(:,1) +shift;

%% Build boundary structure
aux = whos('Tb*');
boundaryNames = {aux.name};
clear aux
refEl = createReferenceElement(1,size(T,2),[]);
nOfBound = numel(boundaryNames);
boundary = struct;
for ibound = 1:nOfBound
    name = boundaryNames{ibound};
    eval(['aux_bound = Tb_' name(4:end) ';'])
    boundary.(name) = aux_bound;
end

Dirichlet_boundaries = setDirichletFaces(boundaryNames);
Nv = size(refEl.NodesCoord,1); 
%% Apply scaling
X = X/lscale;

%% Nodal Connectivity
nNodes = max(max(T));
np = size(T,2);
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
    Te = T(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:np
        N(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) +1;
end
N(:,max(nn):end) = [];

%% Load solution
load('Saves/u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat_Diffusion_0.0001_tau1_Drift_SC2_LocDiffPts2.mat')


figure, plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1),refEl,20,0)
hold on
plotMesh(X*lscale,T)

%% HDG preprocess
disp('HDG preprocess...')
[F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,Dirichlet_boundaries);

%% Shock capturing stuff
shock_st = shockCaptPrep(refEl);
eps_elem = findCoeffShockCaptur(u0,X,T,refEl,0,neq,1e-1);

%% Loop in elements
eps_plot = zeros(size(X,1),1);
xybar = zeros(size(T,1),2);
indel = false(size(T,1),2);
for iElem = 1:size(T,1)
    
    Te = T(iElem,:);
    Xe = X(Te,:);
    % shock capturing parameter in each node
    switch shockcapt
        case 1
            if ~eps_elem(iElem), continue, end
            
            % constant value in each element
            eps_nodal_n = eps_elem(iElem)*ones(Nv,1);
            eps_nodal_u = eps_elem(iElem)*ones(Nv,1);
            eps_plot(Te) = eps_nodal_n;
        case 2
            % linear interpolation
            eps_nodal_n = eps_elem(iElem)*zeros(Nv,1);
            eps_nodal_u = eps_elem(iElem)*zeros(Nv,1);
            els = N(Te(1:3),:);
            aux_eps_nodal_n = zeros(3,1);
            aux_eps_nodal_u = zeros(3,1);
            for i = 1:size(els,1)
                ind = els(i,els(i,:)~=0);
                aux_eps_nodal_n(i) = max(eps_elem(ind));
                aux_eps_nodal_u(i) = max(eps_elem(ind));
            end
            if all(aux_eps_nodal_n==0), continue,end
            eps_nodal_n = shock_st.N*aux_eps_nodal_n;
            eps_nodal_u = shock_st.N*aux_eps_nodal_u;
            eps_plot(Te) = eps_nodal_n;
    end
    xybar(iElem,:) = sum(Xe(1:3,:))/3;
    indel(iElem) = true;
end
figure, plotSolution(X,T,eps_plot,refEl), 
hold on, plotMesh(X*lscale,T);
axis off

%% Circular line
% Density and velocity
rho = u0(1:2:end-1);
gamma = u0(2:2:end);
upar = gamma./rho;
np = 200;
r = 0.75;
theta = pi/2+linspace(0,2*pi,np)';
[line_x, line_y] = pol2cart(theta,r*ones(size(theta)));
% line = ([line_x+shift, line_y] );
line = ([line_x+shift, line_y] )/lscale;
uplot = evalDGapproximationAtPoints(line,upar,X,T,refEl);
figure
plot(theta,uplot,'k-'),grid on
ylabel('Mach')
xlabel('Angle')
xlim([min(theta), max(theta)]);
xybar = xybar(indel,:)*lscale;
xybar(:,1) = xybar(:,1)-shift;
[thetabar,rhobar] = cart2pol(xybar(:,1),xybar(:,2));
thetabar = thetabar-pi/2;
thetabar(thetabar<0) = thetabar(thetabar<0)+2*pi;
thetabar= thetabar+pi/2;
hold on, plot(thetabar,zeros(numel(thetabar)),'ro')
