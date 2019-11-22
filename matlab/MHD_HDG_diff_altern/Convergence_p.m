% HDG method for the Isothermal Euler equations
clear
close all
home

global kT ionmass testcase Mesh diff_n diff_u axisym
kT = 100;
ionmass = 1;


%% Define problem
% problem = 'NS_time';
problem = 'NS_steady';
% problem = 'S_time';
% problem = 'S_steady';

%% Test case
testcase.n = 25;
testcase.wavex = 8;
testcase.wavey = 8;
diff_n = 1;
diff_u = 1;
axisym = 1;

%% Meshes
path2Mesh = '../../Meshes/Meshes_2D/WEST/';
meshName = 'West_h005';

%% convergence points
np = 5; % interpolation degrees

%% number of time steps
nStep = 30;%100000;

%% number of time step to reach the wanted inflow vel
nStep_start = 1;

%% max iteration for Newton Raphson
maxIter = 30;

%% time step
dt0 = 1e-1;

%% base of increment for delta t
baseIncrDt = inf;

%% stability parameter
tau = kT;

%% Iterative solver
iterative = 0;

%% Tolerance Newton-Raphson
tolNR = 1e-6;

%% ************************************************************************
%%
%%                       START COMPUTATIONS
%%
%% ************************************************************************

L2err_dens = zeros(np);
L2err_pvel = zeros(np);
elsize = L2err_dens;
for ip =1:np
    
       %% Prepare case
           load([path2Mesh meshName '_P' num2str(ip) '.mat'])
    aux = whos('Tb*');
    boundaryNames = {aux.name};
    refEl = createReferenceElement(1,size(T,2),[]);
    nOfBound = numel(boundaryNames);
    boundary = struct;
    for ibound = 1:nOfBound
        name = boundaryNames{ibound};
        eval(['aux_bound = Tb_' name(4:end) ';'])
        boundary.(name) = aux_bound;
    end
    

        %% Apply translation if axisymmetric case
%         if (axisym && min(X(:,1))<0) || testcase.n==25 
%                 X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+1;
%         end

        % Fill global variable Mesh
        Mesh.X = X;
        Mesh.T = T;
        Mesh.lscale = 1;
        for ibound = 1:nOfBound
            name = boundaryNames{ibound};
            eval(['aux_bound = Tb_' name(4:end) ';'])
            Mesh.(name) = aux_bound;
        end
        Mesh.maxx = max(X(:,1));
        Mesh.minx = min(X(:,1));
        Mesh.maxy = max(X(:,2));
        Mesh.miny = min(X(:,2));
        
        
        % set which faces are Dirichlet and the pressure condition
        Dirichlet_boundaries = setDirichletFaces(boundaryNames);
        
        % HDG preprocess
        disp('HDG preprocess...')
        [F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,Dirichlet_boundaries);
        
        % dimensions
        neq = 2;                         % n of equations (2 for MHD isothermal 2D)
        nf = size(F,2);                  % n of faces per element
        Nv = size(refEl.NodesCoord,1);   % n of element nodes for the velocity
        nv = size(refEl.faceNodes1d,2);  % n of face nodes for the velocity
        Nf_dir = sum(sum(F_dir));         % n of faces on the Dirichlet boundary
        Nf = max(max(F));                 % n of faces
        Nf_unk = Nf-Nf_dir;               % n of faces without the Dirichlet boundary
        nvUnk = neq*nv*Nf_unk;              % velocity unknowns
        nvTot = neq*nv*Nf;                  % velocity nodes on the skeleton
        Ne = size(T,1);                   % n of elements
        nIntFaces = size(infoFaces.interiorFaces,1); % number of interior faces
        
        % element size
%         elsize(ic,ip) = 1/2^ic;
        
        % precalculated matrices
        [M,B,C,C_dir,G,D,E,Df,Ef,Edir,invL,force,Lf,...
            L,P,Pb,Q,Qb,Qf] =...
            hdg_PrecalculatedMatrices(X,T,flipFace,refEl,tau,F,F_dir);
        
        % initialize convection matrices
        Cv = zeros(neq*Nv,neq*Nv,Ne);
        H = zeros(neq*Nv,nf*neq*nv,Ne);
        Hf = zeros(nf*neq*nv,nf*neq*nv,Ne);
        Hdir = zeros(neq*Nv,Ne);
        Hfdir = zeros(nf*neq*nv,Ne);
        fH = Hfdir;
        % check the problem type
        checkProblem
        
        % initialization
        % u0 = initializeSolution(X,T,F_dir,refEl);
        u0 = initializeSolutionToAnalyticalSol(X,T);
        u_tilde = extractFaceSolution(u0,F,F_dir,refEl);
        u_tilde = [u_tilde; zeros(nvTot-nvUnk,1)];
        dt = dt0;
        time = 0;
        tic
        %% Time integration
        for iStep = 1:nStep
            disp(['Step: ' num2str(iStep)])
            
            time = time + dt;
            u_iter = u0;
            
            %% Newton Raphson iterations
            for iter = 1:maxIter
                
                % Calculate convection matrices (only if NS)
                if strcmp(problem,'NS_time') || strcmp(problem,'NS_steady')
                    [Cv,H,Hdir,Hf] = ...
                        hdg_ConvectionMatrices(X,T,F,flipFace,refEl,u_iter,u_tilde,F_dir);
                end
                
                % Get the mapping
                [LL,LL0,UU,UU0] = hdg_Mapping(flipFace,Nv,nv,M,G,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
                    B,C,C_dir,invL,P,Pb,Q,Qb);
                
                % Impose boundary condition of the global problem
               [Df,Ef,Hf,Lf,Qf,fH] = hdg_bc(Df,Ef,Hf,Lf,Qf,X,T,infoFaces,refEl,boundaryNames,F,u_tilde);
                
               % Assembly
               [K,f] = hdg_Assembly(F,nv,LL,LL0,UU,UU0,Df,Ef,Lf,...
                   Hf,Qf,fH,F_dir,nIntFaces);
                
                % Face solution
                sol = K\f;
                u_tilde = [sol(1:nvUnk); zeros(nvTot-nvUnk,1)];
                
                % element by element solution
                u = calculateElementByElementSolution(u_tilde,F,UU,UU0);
                
                
                % check convergence for Newton Raphson
                errorNR = norm(u-u_iter);
                if errorNR<tolNR
                    disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
                    break
                elseif errorNR > 1e4
                    error('Problem in the N-R procedure')
                else
                    u_iter = u;
                    disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
                end
                
            end
            
            disp(['Newton Raphson iterations: ' num2str(iter)])
            if iter == maxIter && iter~=1
                error('Problem: Newton Raphson does not converge')
            end
            
            % check convergence for the time iterations
            errorTime = norm(u-u0);
            if errorTime < 1e-4
                break
            else
                u0 = u;
                dt = dt0*2^(-log(errorTime)/log(baseIncrDt));
                disp(['Error in time step '  num2str(iStep) ': ' num2str(errorTime)])
                disp(['Time step: ' num2str(dt)])
                disp(['Simulation time: ' num2str(time)])
                disp(' ')
            end
        end
    
    L2err = computeErrorAnalyticSol(X,T,u,refEl);
    L2err_dens(ip) = L2err(1);
    L2err_pvel(ip) = L2err(2);
end

save(['Results/Conv_p_Np' num2str(np) ...
    'Testcase' num2str(testcase.n) '_Wx' num2str(testcase.wavex) ...
    '_Wy' num2str(testcase.wavey) '.mat'], 'L2err_dens', 'L2err_pvel')
