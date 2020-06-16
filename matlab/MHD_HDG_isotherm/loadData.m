%% Load meshes 
% X = nodal coordinates for velocity
% T = connectivity for velocity
% Tb* boundary connectivity for velocity
% 
% %  
meshName = 'mesh2_P4.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

%  meshName = 'Circle_lim_h02_P4.mat';
% % meshName = 'circ_2_P2.mat';
% path2Mesh = '../../Meshes/Meshes_2D/';


% meshName = 'West_NoHole_h0.04_refCorn_P4.mat';
% meshName = 'West_SmallHole_h0.04_refCorn_P3.mat';
% meshName = 'West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P2.mat';
% meshName = 'West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P2.mat';
% path2Mesh = '../../Meshes/Meshes_2D/WEST/';
 
% meshName = 'Circle_Hole_New_h0.1_P4.mat';
% path2Mesh = '../../Meshes/';
% % 
% meshName = 'Circle_lim_h01_refCor_P2.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.006_refSep0.008_P6.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.83_refCorn0.001_refSol0.004_refSep0.008_P3.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P5.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.006_refSep0.002_DoubleSol_P8.mat';
% meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
% meshName = 'Circel_LIM_small_RefCircle_h0.01_h0.1_P2.mat';
% meshName = 'Circle_lim_h01_P8.mat';
% meshName = 'Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P2.mat';

%  meshName = 'Circle_CAMILLE_P2.mat';
% path2Mesh = '../../Meshes/Meshes_2D/';


% meshName = 'Cylinder_P2.mat';
% path2Mesh = '../../Meshes/Meshes_NS/';

%% Lenght scale
if testcase.n>=50 && testcase.n<100
    lscale = 1.901*1e-3;
%     lscale = 2;
else
    lscale = 1;
end

%% load Mesh
if strcmpi(meshName(end-2:end),'dcm')
    GenerateMatFileFromEZ4U([path2Mesh meshName]);
    load([path2Mesh meshName(1:end-4) '_readed.mat'])
    delete([path2Mesh meshName(1:end-4) '_readed.mat'])
elseif strcmpi(meshName(end-2:end),'mat')
    load([path2Mesh meshName])
else
    error('wrong name')
end

%% Set Nemann faces (if needed)
if setneum
   setNeumannFaces; 
end

%% Apply translation if axisymmetric case
if (axisym && min(X(:,1))<1e-8) || testcase.n==25
    
    if testcase.n==25
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+0.6;
    else
        % apply translation in x
%         X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+3.4;
         X(:,1) = X(:,1) +3.4;
    end
end

%% Apply scaling
X = X/lscale;


if dirichlet_weak
   Tb_Dweak = Tb_Diriclet; 
   clear Tb_Diriclet
   elementFaceInfo.Dweak = elementFaceInfo.Diriclet;
   elementFaceInfo = rmfield(elementFaceInfo,'Diriclet');
end
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

%% Fill global variable Mesh
Mesh.X = X;
Mesh.T = T;
Mesh.boundaryNames = boundaryNames;
for ibound = 1:nOfBound
    name = boundaryNames{ibound};
    eval(['aux_bound = Tb_' name(4:end) ';'])
    Mesh.(name) = aux_bound;
end
Mesh.lscale = lscale;
Mesh.maxx = max(X(:,1));
Mesh.minx = min(X(:,1));
Mesh.maxy = max(X(:,2));
Mesh.miny = min(X(:,2));


% Nodal Connectivity
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
Mesh.N = N;
