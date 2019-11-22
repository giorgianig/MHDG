% HDG method for the Isothermal Euler equations
clear
clear global
close all
home

global kT testcase Mesh refEl refElTor neq diff_n diff_u Magnetic coreflux axisym
global driftvel  shockcapt localDiffPoints limitRhoMin useThreshold
global ntor ptor theta dirichlet_weak coeffconv
kT = 25;%12;
% kT =  0.024;
neq = 2;   % n of equations (2 for MHD isothermal 2D)



axisym = 1;
coeffconv = 1;


driftvel = 0;
shockcapt = 0; % 0- no sh 1-elementwise 2-linear interp
localDiffPoints = 0;
limitRhoMin = 0;
useThreshold = 0;
MultAddDiff = 1;
dirichlet_weak = 1;

%% Define problem
% problem = 'NS_time';
problem = 'NS_steady';
% problem = 'NS_pseudotime';
% problem = 'S_steady';

%% Test case
testcase.n = 2;
testcase.wavex = 1;
testcase.wavey = 0;
testcase.xc = -0.5;
testcase.yc = -0.5;
setneum = 0; % Neumann BC
diff_n = 3;%3.8;
diff_u = 3;%3.8;
coreflux = [1;0];
ntor = 2;                % number of toroidal division
ptor = 1;                 % degree in the poloidal direction
theta = 2*pi;         % Toroidal angle considered
ntpos =3;






%% load data
loadData

%% number of step for each plot
nPlot = 1;

%% number of time steps
nStep = 100000;

%% number of time step to reach the wanted inflow vel
nStep_start = 1;

%% max iteration for Newton Raphson
maxIter = 4;

%% time step
dt0 = 1e9;

%% base of increment for delta t
baseIncrDt = inf;

%% stability parameter
% tau = kT/Mesh.lscale;
tau = 1;

%% Iterative solver
iterative = 0;

%% Tolerance Newton-Raphson
tolNR = 1e-6;

%% Tolerance for the time convergence
tolTime = 1e-8;

%% Restart
restart = 0;

%% Pseudo time step to lower diffusion
diffscale = 1/10^(0.5); % diffusion is multiplied by this value in each time step
% in case of NS_pseudotime
minDiff = 1e-6;  % minimum value of diffusion in case pseudotime
psStepChange = 1;
%% ************************************************************************
%%
%%                       START COMPUTATIONS
%%
%% ************************************************************************

dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);


refElTor = createReferenceElement(0,(ptor+1)^2,[]);
% Creating faceNodes3 for face numbering in the prismatic elements
Np2d = size(refEl.NodesCoord,1);           % Number of points for 2d element
faceNodes3 = zeros(size(refEl.faceNodes,1),size(refEl.faceNodes,2)*size(refElTor.faceNodes1d,2));
for ifa = 1:size(refEl.faceNodes,1)
    faceNodes3(ifa,:) = transpose(col(bsxfun(@plus,refEl.faceNodes(ifa,:)' ,(refElTor.faceNodes1d-1)*Np2d   )));
end

% Creating basis functions in the toroidal faces for prismatic elements
Np1dPol = size(refEl.N1d,2);
Np1dTor = size(refElTor.N1d,2);
Ng1dPol = size(refEl.N1d,1);
Ng1dTor = size(refElTor.N1d,1);
shapeFunctionsTF = zeros(Ng1dPol*Ng1dTor,Np1dPol*Np1dTor);
for i = 1:Np1dTor
    for j=1:Np1dPol
        k = (i-1)*Np1dPol+j;
        shapeFunctionsTF(:,k) = col(refEl.N1d(:,j)*refElTor.N1d(:,i)');
    end
end

% stop
refElTor.sfTor = shapeFunctionsTF;
refElTor.faceNodes3 = faceNodes3;





%% Initialization
if restart
    %     load('../MHD_HDG_diff_altern/Saves/u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_0.0031623_tau1_UseTh1e-06.mat','u0','u_tilde');
    load('../MHD_HDG_diff_altern/Saves/u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat_Diffusion_0.1_tau1.mat','u0','u_tilde');
    
    %         u0 = projectFieldDifferentP(u0,4,refEl.degree,size(T,1));
else
    % u0 = initializeSolution(X,T,F_dir,refEl);
    u0 = initializeSolutionToAnalyticalSol_3D(X,T,refEl);
    % u0 = initializeSolutionWithL2Projection(X,T,refEl);
end
%
% figure, plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,5)
% stop

% u0 = projectSolutionDifferentMeshes('Circle_lim_h01_refCor_P10.mat','u_LIM_h01_refCorn_diff0.02_P10.mat');






%% set which faces are Dirichlet
Dirichlet_boundaries = setDirichletFaces(boundaryNames);
% Dirichlet_boundaries = false(size(Dirichlet_boundaries));


%% HDG preprocess
disp('HDG preprocess...')
[F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,Dirichlet_boundaries);

if exist('u_tilde','var')
    %         u_tilde = projectUtildeDifferentP(u_tilde,4,refEl.degree,max(max(F)));
end

% Load magnetic field if needed
if (testcase.n >= 50 && testcase.n<60)
    disp('Loading magnetic field...')
    loadMagneticField(X,T,F,refEl);
end



%% Set dimensions
nf = size(F,2);                  % n of faces per element
Nf_dir = sum(sum(F_dir));         % n of faces on the Dirichlet boundary
Nf = max(max(F));                 % n of faces
Nf_unk = Nf-Nf_dir;               % n o...f faces without the Dirichlet boundary
nIntFaces = size(infoFaces.interiorFaces,1); % number of interior faces
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nfl   = Np1dPol*Np1dTor;
Np2D  = size(refEl.NodesCoord,1);
N2D   = size(T,1);
nvTot  = ntor*neq*(Nfl*Nf + Np2D*N2D);
nvUnk  = ntor*neq*(Nfl*Nf_unk + Np2D*N2D);
% element size
h = 1;

%% Precalculated matrices
disp('Precalculated matrices...')
[M,B,C,C_dir,G,D,E,Df,Ef,Edir,invL,force,Lf,...
    L,P,Pb,Q,Qb,Qf] =...
    hdg_PrecalculatedMatrices3D(X,T,flipFace,refEl,tau,F,F_dir);



%% Uncomment these lines for isotropic diffusion
% Qf = zeros(size(Qf));
% Pb = zeros(size(Pb));
% Qb = zeros(size(Qb));

%% Initialize convection matrices
% Cv = zeros(neq*Nv,neq*Nv,Ne);
% H = zeros(neq*Nv,nf*neq*nv,Ne);
% Hf = zeros(nf*neq*nv,nf*neq*nv,Ne);
% Hdir = zeros(neq*Nv,Ne);
% Hfdir = zeros(nf*neq*nv,Ne);
% fH = Hfdir;


%% check the problem type
checkProblem

%% Extract u tilde
if ~exist('u_tilde','var')
    u_tilde = extractFaceSolution_3d(u0,F,F_dir,infoFaces,refEl);
    u_tilde = [u_tilde; zeros(nvTot-nvUnk,1)];
end

%% Some graphics before starting the computation
% figure, plotSolution(X,T,u0(1:2:end-1),refEl,20); title('U(1)')
% figure, plotSolution(X,T,u0(2:2:end),refEl,20); title('U(2)')
% GradPerp = computeGradPerp(X);
% GradPerp = transpose(reshape(GradPerp,[4,numel(GradPerp)/4]));
% figure, plotSolution(X,T,GradPerp(:,1),refEl,20); title('dU(1)/dx')
% figure, plotSolution(X,T,GradPerp(:,2),refEl,20); title('dU(1)/dy')
% figure, plotSolution(X,T,GradPerp(:,3),refEl,20); title('dU(2)/dx')
% figure, plotSolution(X,T,GradPerp(:,4),refEl,20); title('dU(2)/dy')
% stop
% subplot(1,2,1),plotSolution(X,T,u0(1:2:end-1),refEl); title('Density HDG')
% subplot(1,2,2),plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
% drawnow

store_simulation_parameters;

%% Time integration
dt = dt0;
time = 0;
psStep = 0; % pseudostep
tic
for iStep = 1:nStep
    disp(['Step: ' num2str(iStep)])
    
    time = time + dt;
    u_iter = u0;
    psStep = psStep+1;
    
    ifig=0;
    %% Newton Raphson iterations
    for iter = 1:maxIter
        
        ifig = ifig+1;
        
        % Calculate convection matrices (only if NS)
        if strcmp(problem,'NS_time') || strcmp(problem,'NS_steady') || strcmp(problem,'NS_pseudotime')
            [Cv,H,Hdir,Hf] = ...
                hdg_ConvectionMatrices_3D(X,T,F,flipFace,refEl,u_iter,u_tilde,F_dir);
        end...
        
    if coeffconv==0
        Cv = Cv*0;
        H = H*0;
        Hdir = Hdir*0;
        Hf = Hf*0;
        G = G*0;
    end
    
%         residual = computeResiduals_local_analiticalGrad_3d...
%     (F,F_dir,flipFace,u_tilde,u0,M,G,Cv,D,E,Edir,H,Hdir,force,...
%     B,C,C_dir,L,P,Pb,Q,Qb,invL,refEl,infoFaces)
% stop

        % Get the mapping
        [LL,LL0,UU,UU0] = hdg_Mapping3d(flipFace,refEl,refElTor,M,G,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
            B,C,C_dir,invL,P,Pb,Q,Qb);
                
%         checkMapping3d(LL,LL0,UU,UU0,u_tilde,F,F_dir,T,X)

% LLfort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/LL.txt');
% UUfort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/UU.txt');
% L0fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/L0.txt');
% U0fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/U0.txt');
% 
% 
% max(abs(LLfort(:)-LL(:)))
% max(abs(UUfort(:)-UU(:)))
% max(abs(L0fort(:)-LL0(:)))
% max(abs(U0fort(:)-UU0(:)))
        % Impose boundary condition of the global problem
        [Df,Ef,Hf,Lf,Qf,fH] = hdg_bc3d(Df,Ef,Hf,Lf,Qf,X,T,infoFaces,refEl,boundaryNames,F,F_dir,u_tilde);
        
        %          Hf = zeros(size(Hf));
% disp('Checking')        
 debugFortran
                       stop
        %                stop
        % checkFortMapping
        % %
        % stop
        % Assembly
        [K,f] = hdg_Assembly3d(refEl,F,F_dir,LL,LL0,UU,UU0,Df,Ef,Lf,Hf,Qf,fH);
        
        %  matK_ser = readCSRtxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_parall/MatK_0_1.txt');
        %  matK = readCSRtxtParall('/home/giorgio/Dropbox/Fortran/MHDG/test/test_parall/MatK',2);
        %  f_fort = readCSRVectortxtParall('/home/giorgio/Dropbox/Fortran/MHDG/test/test_parall/rhs',2);
        % % % f_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/fRHS.txt');
        % figure, spy(matK)
        % figure, spy(K)
        % disp(['Error in global matrix: ' num2str(max(max(abs(matK-K))))])
        % % % disp(['Error in global RHS: ' num2str((max(abs(f-f_fort))))])
        % stop
        
        %         Hf = zeros(size(Hf));
        %         [Kaux,faux] = hdg_Assembly(F,nv,U,U0,Df,Ef,Hf,fH,F_dir,nIntFaces);
        %         max(max(abs(K-Kaux)))
        %         stop
        %         condest(K)/1e9
        %         stop
        % face solution
        sol = K\f;
        
        
        
        u_tilde = [sol(1:nvUnk); zeros(nvTot-nvUnk,1)];
        
        % element by element solution
        u       = computeElementByElementSolution3d(ntor,u_tilde,UU,UU0,T,F,F_dir,refEl,refElTor);
        gradU   = computeElementByElementGradient3d(ntor,u_tilde,LL,LL0,T,F,F_dir,refEl,refElTor);
        
        % figure(ifig)
        % subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
        % subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
        % drawnow
        %% Checking




        %                 stop
%         resglob  = computeResiduals_global...
%             (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Hf,fH,F_dir);
%         sum(sum(resglob))
%         resglob_alt  = computeResiduals_global_alt...
%             (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Hf,fH,...
%             LL,LL0,UU,UU0,F_dir);
%         
%         checkMapping(LL,LL0,UU,UU0,u_tilde,F,T,X);
%         residual = computeResiduals_local_3d...
%             (F,F_dir,flipFace,u_tilde,u,gradU,M,G,Cv,D,E,Edir,H,Hdir,force,...
%     B,C,C_dir,L,P,Pb,Q,Qb,invL,refEl,infoFaces)
%         sum(sum(residual))
%         stop
        
        %% Check and fix problems in the solution
        if useThreshold~=0
            u = checksolution(u);
        end
        
        % check convergence for Newton Raphson
        errorNR = norm(u-u_iter)/sqrt(numel(u));
        if errorNR<tolNR
            disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
            break
        elseif errorNR > 1e6
            error('Problem in the N-R procedure')
        else
            u_iter = u;
            disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
        end
        
    end
    
    disp(['Newton Raphson iterations: ' num2str(iter)])
    %     if iter == maxIter && iter~=1
    %         error('Problem: Newton Raphson does not converge')
    %     end
    
    % plot
    if ~mod(iStep,nPlot)
        converged_case = 0;
        saveSolutionWithParametersName;
        %         plotDiscontinuousSolution(X,T,u(1:3:end-2),refEl);
        % subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl,20,0,1); title('Density HDG')
        % subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
        
        %         storeSolution(X,T,Xp,Tp,L_gradv,u,p,refEl,boundary)
        %                 drawnow
    end
    
    
    %     stop
    
    
    % check convergence for the time iterations
    errorTime = norm(u-u0)/sqrt(numel(u));
    if errorTime < tolTime
        if strcmpi(problem,'NS_time')
            disp(['Error in time ' num2str(errorTime)])
            break
        elseif strcmpi(problem,'NS_pseudotime')
            disp(['Error in time ' num2str(errorTime)])
            if diff_n>minDiff && diff_u>minDiff
                if psStep>psStepChange % We make at least 10 time step with this diffusion
                    disp(['*********Changing diffusion**************  '  num2str(diff_n*diffscale)])
                    u0 = u;
                    store_simulation_parameters;
                    converged_case = 1;
                    saveSolutionWithParametersName;
                    Lf = Lf*diffscale;
                    Qf = Qf*diffscale;
                    P  = P*diffscale;
                    Pb = Pb*diffscale;
                    Q  = Q*diffscale;
                    Qb = Qb*diffscale;
                    diff_n = diff_n*diffscale;
                    diff_u = diff_u*diffscale;
                    psStep = 0;
                    dt = dt0;
                else
                    u0 = u;
                    dt = dt0*2^(-log(errorTime)/log(baseIncrDt));
                    disp(['psStep: ' num2str(psStep)])
                    disp(['Time step: ' num2str(dt)])
                end
            else
                store_simulation_parameters;
                saveSolutionWithParametersName;
                break
            end
        end
    else
        u0 = u;
        dt = dt0*2^(-log(errorTime)/log(baseIncrDt));
        disp(['Error in time step '  num2str(iStep) ': ' num2str(errorTime)])
        disp(['Time step: ' num2str(dt)])
        disp(['Simulation time: ' num2str(time)])
        disp(' ')
    end
    
    
end
toc

store_simulation_parameters;
converged_case = 1;
saveSolutionWithParametersName;

ur = transpose(reshape(u,[neq,numel(u)/neq]));
% Plot solution
disp('Plotting...')
for itor = 1:ntpos
    %     upol = extractSolutionInAPoloidalPlane(u,T,refElPol,refElTor,itor);
    upol = extractSolutionInAtGivenTheta(ur(:,1),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - density'])
    upol = extractSolutionInAtGivenTheta(ur(:,2),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - momentum'])    
end


%% Plot analytical solution
disp('Plotting analytical solution...')
th = linspace(0,theta,ntor+1);
uan = initializeSolutionToAnalyticalSol_3D(X,T,refEl);
uan = transpose(reshape(uan,[neq,numel(uan)/neq]));
for itor = 1:ntpos
    %      upol = analyticalSolution(X,th(itor));
    upol = extractSolutionInAtGivenTheta(uan(:,1),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Anal sol plane ' num2str(itor) ' - density'])
    upol = extractSolutionInAtGivenTheta(uan(:,2),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Anal sol plane ' num2str(itor) ' - momentum'])    
end



%% L2 error
disp('Computing error...')
err = computeError_3D(u,X,T,theta,refEl,refElTor,ntor);
disp(['Error U1 = ', num2str(err(1),10)]);
disp(['Error U2 = ', num2str(err(2),10)]);
disp(' ')

