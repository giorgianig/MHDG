clear all
close all

meshData = 'square_h0.1';

% generate mesh
GenerateMatFileFromEZ4U(['Meshes/' meshData '_P8.dcm']);
load(['Meshes/' meshData '_P8_readed.mat'])
delete(['Meshes/' meshData '_P8_readed.mat'])

% plotMesh
plotMesh(X,T);

% analytical solution
u = analyticalSolution(X);

% velocity plot
hold on 
quiver(X(:,1),X(:,2),u(:,1),u(:,2));
axis equal;
