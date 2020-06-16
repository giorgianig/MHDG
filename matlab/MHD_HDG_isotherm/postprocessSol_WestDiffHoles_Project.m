% Postprocess routine

clear
close all

global testcase Mesh axisym neq

%% Physical parameters
lscale = 1.901*1e-3;
kT = 1;
neq = 2;
testcase.n = 50;
axisym = 1;
fref = 7.28e6; % reference frequency
D = 0.02;
nref = 5;
npoints = 100;
%% Find the separatrix
load ../../Meshes/Meshes_2D/WEST/WEST_wall.mat
load ../../Meshes/Meshes_2D/WEST/WEST_far_465.mat

%% Load mesh
path2Mesh = '../../Meshes/Meshes_2D/WEST/';
path2Sol = 'Saves/';

%%*************************************************************************
%%
%% Large hole
%%
%%*************************************************************************
meshName = 'West_h0.04_refCorn_P8_fixed.mat';
path2Sol = './Saves/';
load([path2Mesh meshName]);

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

% CASE WITHOUT DRIFT
% solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51.mat';
% load([path2Sol solName])
% % figure,plotSolution(X,T,u0(1:2:end-1),refEl,nref)
% % Project the solution without hole on the mesh with large hole
%  meshNoHole = [path2Mesh 'West_NoHole_h0.04_refCorn_P8.mat'];
%  solNoHole = [path2Sol 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51.mat'];
% figure, plotSolutionDifference(X,T,u0(1:2:end-1),refEl,solNoHole,meshNoHole,nref)
% 
% % upr = projectSolutionDifferentMeshes(meshNoHole,solNoHole);
% % upr = col(upr');
% % save('Saves/u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_PROJECTEDONLARGEHOLE.mat','upr')
% drawnow


% CASE WITH DRIFT
solName = 'u_West_h0.04_refCorn_P8_fixed.mat_Diffusion_0.02_tau1_tc51_Drift.mat';
load([path2Sol solName])
meshNoHole = [path2Mesh 'West_NoHole_h0.04_refCorn_P8.mat'];
solNoHole = [path2Sol 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat'];
figure, plotSolutionDifference(X,T,u0,refEl,2,solNoHole,meshNoHole,nref)

% Project the solution without hole on the mesh with large hole
% meshNoHole = [path2Mesh 'West_NoHole_h0.04_refCorn_P8.mat'];
% solNoHole = [path2Sol 'u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift.mat'];
% upr = projectSolutionDifferentMeshes(meshNoHole,solNoHole);
% upr = col(upr');
% save('Saves/u_West_NoHole_h0.04_refCorn_P8.mat_Diffusion_0.02_tau1_tc51_Drift_PROJECTEDONLARGEHOLE.mat','upr')


