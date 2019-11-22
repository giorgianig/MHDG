% HDG method for the N-Gamma-Ti-Te model for plasma-edge in 3D
clear all
% clear global
close all
home

global  testcase Mesh refEl neq diff_n diff_u diff_ei diff_ee diff_pari diff_pare Magnetic coreflux axisym epn tau
global driftvel  shockcapt localDiffPoints limitRhoMin useThreshold Gmbohmi Gmbohme Mref tie pcourr prolongateExponents
global useNeumannBC decoupleEquations physValueGar decoupleTiTe constantExchangeTemp
global refElTor 
global ntor ptor theta dirichlet_weak


neq = 4;   % n of equations (4 N-Gamma-Ti-Te model)

driftvel = 0;
shockcapt = 0; % 0- no sh 1-elementwise 2-linear interp
localDiffPoints =0;
limitRhoMin = 0;
useThreshold =0;
MultAddDiff = 1;
dirichlet_weak = 1;

%% Define problem
% problem = 'NS_time';
problem = 'NS_steady';
% problem = 'NS_pseudotime';
% problem = 'S_steady';

%% Test case
testcase.n =60;
axisym = 1;
testcase.wavex = 1;
testcase.wavey = 1;
testcase.xc = inf;
testcase.yc = inf;
setneum = 0; % Neumann BC
diff_n = 0.38;
diff_u = 0.38;
diff_ei = 0.38;
diff_ee = 0.38;
diff_pari = 3e5;%20;
diff_pare = 1e7;%30;
epn = 5/2;
tie = 1e10;
coreflux = [1;0; 18;18];
Gmbohmi = 2.5;
Gmbohme = 4.5;
Mref = 12;
ntor = 1;                % number of toroidal division
ptor = 1;                 % degree in the poloidal direction
theta = 6.28;         % Toroidal angle considered
ntpos =1;
dumpNR = 0.25;


pcourr = 1; % only for debug: should be 1
prolongateExponents = 0;
useNeumannBC = 0;        % use a Neumann BC instead of Bohm (but still sonic speed)
% 0- Bohm classic
% 1- Grad U = 0 + sonic speed
% 2- Grad// T = 0 + sonic speed
decoupleEquations =0;    % change Jacobian in order to decouple n-Gamma from Ti-Te
physValueGar =0;
decoupleTiTe =0;            % decouple the two temperature equations
constantExchangeTemp =0;

savememory = 0;

drawsol = 0; % draw sol at NR iterations

%% load data
loadData

%% number of step for each plot
nPlot = 10;

%% number of time steps
nStep = 100000;

%% number of time step to reach the wanted inflow vel
nStep_start = 1;

%% max iteration for Newton Raphson
maxIter = 3;

%% time step
dt0 = 1e-1;

%% base of increment for delta t
baseIncrDt = inf;

%% stability parameter
% tau = kT/Mesh.lscale;
tau = 100;%diff_pare;

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



% figure, plotMesh(X*Mesh.lscale,T)
% stop
% some consistency check
if useNeumannBC && testcase.n<49
    error('Attention, are you sure you want Neumann bc?')
end
% if useNeumannBC && tau~=1
%     error('Attention, you probably want to use tau=1')
% end

%% Initialization
if restart
%         load('../MHD_HDG_diff_altern/Saves/u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_0.0031623_tau1_UseTh1e-06.mat','u0','u_tilde');
        load('Saves/u_West_Hugo_h0.04_refCorn_P3.mat_Dpe10_Dpai1000_Dpae1000_Mref12_tau20_tc50_UseTh0.1.mat','u0','u_tilde','gradU');
    %         load('Saves/u_Circle_limInfThin_h0.1_P8.mat_Dn10_Dg10_De10_Dpar0_Mref12_tau1000_tc60_UseTh1e-06.mat','u0','u_tilde','gradU');
    %         load('Saves/Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P4_Dper10_Dpar0_Neum.mat','u','u_tilde','gradU'); u0=u;
    
%     HDF5load(['/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/'...
%         'Sol_Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P2_Diff.20000E+01_CONV.h5'])
    
%     gradU = q';
%     u0 = u';
%     u_tilde = u_tilde';
    
    % clear q
    %     gradU = reshape(gradU,[2,numel(gradU)/2]);
    %     u0x = gradU(:,1);
    %     u0y = gradU(:,2);
    %         u0 = projectFieldDifferentP(u0,4,refEl.degree,size(T,1));
else
    % u0 = initializeSolution(X,T,F_dir,refEl);
    [u0,u0x,u0y,u0t] = initializeSolutionToAnalyticalSol_3D(X,T,refEl);
    gradU = col(transpose([u0x,u0y,u0t]));
    % u0 = initializeSolutionWithL2Projection(X,T,refEl);
end
%

% solq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/solq.txt');
% max(abs(gradU-solq))
% solu = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/solu.txt');
% max(abs(u0-solu))
% solutilde = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/solutilde.txt');
% max(abs(u_tilde-solutilde))
% stop


% figure, plotSolution(X,T,u0(2:2:end)./u0(1:2:end-1)/sqrt(kT),refEl,5)
% stop

% u0 = projectSolutionDifferentMeshes('Circle_lim_h01_refCor_P10.mat','u_LIM_h01_refCorn_diff0.02_P10.mat');

up = cons2phys(transpose(reshape(u0,[neq,numel(u0)/neq])));
ifig = 0;

% stop


%% set which faces are Dirichlet
Dirichlet_boundaries = setDirichletFaces(boundaryNames);
% Dirichlet_boundaries = false(size(Dirichlet_boundaries));


%% HDG preprocess
disp('HDG preprocess...')
[F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,refEl,Dirichlet_boundaries);

if exist('u_tilde','var')
    %         u_tilde = projectUtildeDifferentP(u_tilde,4,refEl.degree,max(max(F)));
end

% Load magnetic field if needed
if (testcase.n >= 50 && testcase.n<60)
    disp('Loading magnetic field...')
    loadMagneticField(X,T,F,refEl);
end

if shockcapt
    shock_st = shockCaptPrep(refEl);
end


%% Set dimensions
%% Set dimensions
nf = size(F,2);                  % n of faces per element
Nf_dir = sum(sum(F_dir));         % n of faces on the Dirichlet boundary
Nf = max(max(F));                 % n of faces
Nf_unk = Nf-Nf_dir;               % n of faces without the Dirichlet boundary
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
[M,B,C,C_dir,D,E,Df,Ef,Edir,invL,force,Lf,...
    L,P,Pb,Q,Qb,Qf] =...
    hdg_PrecalculatedMatrices3D(X,T,flipFace,refEl,tau,F,F_dir);

if localDiffPoints
    [Lf_locDiffPts,P_locDiffPts,Q_locDiffPts] = hdg_LocalDiffusionInPoints(X,T,flipFace,refEl,Lf,P,Q);
else
    Lf_locDiffPts = 0;
    P_locDiffPts = 0;
    Q_locDiffPts = 0;
end

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
if drawsol
    figure(100+ifig)
    subplot(2,5,1),plotSolution(X,T,up(:,1),refEl,5);  title('Density')
    subplot(2,5,2), plotSolution(X,T,up(:,2),refEl,5); title('Parallel velocity')
    subplot(2,5,3),plotSolution(X,T,up(:,3),refEl,5);  title('Total energy for ions ')
    subplot(2,5,4), plotSolution(X,T,up(:,4),refEl,5); title('Total energy for electrons')
    subplot(2,5,5), plotSolution(X,T,up(:,5),refEl,5); title('Ions pressure HDG')
    subplot(2,5,6), plotSolution(X,T,up(:,6),refEl,5); title('Electrons pressure HDG')
    subplot(2,5,7), plotSolution(X,T,up(:,7),refEl,5); title('Ions temperature HDG')
    subplot(2,5,8), plotSolution(X,T,up(:,8),refEl,5); title('Electrons temperature HDG')
    subplot(2,5,9), plotSolution(X,T,up(:,9),refEl,5); title('Sound speed HDG')
    subplot(2,5,10), plotSolution(X,T,up(:,10),refEl,5); title('Mach')
    % % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    drawnow
    
    
%     subplot(2,2,1),plotSolution(X,T,u0(1:4:end-3),refEl); title('U1')
%     subplot(2,2,2),plotSolution(X,T,u0(2:4:end-2),refEl); title('U2')
%     subplot(2,2,3),plotSolution(X,T,u0(3:4:end-1),refEl); title('U3')
%     subplot(2,2,4),plotSolution(X,T,u0(4:4:end),refEl);    title('U4')
%     drawnow
    
end

% stop
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
    if shockcapt
        [Lf_sc,P_sc,Q_sc] = hdg_ShockCapturing(X,T,u_iter,flipFace,refEl,shock_st);
    else
        Lf_sc = 0;
        P_sc = 0;
        Q_sc = 0;
    end
    if limitRhoMin~=0
        [Lf_limRho,P_limRho,Q_limRho] = hdg_AddDiffForLimitingRho...
            (X,T,F,u_iter,u_tilde,flipFace,refEl,MultAddDiff);
    else
        Lf_limRho = 0;
        P_limRho = 0;
        Q_limRho = 0;
    end
    
    %% Newton Raphson iterations
    for iter = 1:maxIter
        
        ifig = ifig+1;
        
        % Calculate convection matrices (only if NS)
        if strcmp(problem,'NS_time') || strcmp(problem,'NS_steady') || strcmp(problem,'NS_pseudotime')
            [Cv,H,Hdir,Hf] = ...
                hdg_ConvectionMatrices_3D(X,T,F,flipFace,refEl,u_iter,u_tilde,F_dir);

            [ TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf,TUhdir] =...
                hdg_ParallDiffusionEnergyMatrix_3D(X,T,F,flipFace,refEl,gradU,u_iter,u_tilde,F_dir);
            
        end
        
        

        % Get the mapping
        [LL,LL0,UU,UU0] = hdg_Mapping_3D(flipFace,refEl,refElTor,M,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
            B,C,C_dir,invL,P+P_locDiffPts+P_sc+P_limRho,Pb,Q+Q_locDiffPts+Q_sc+Q_limRho,Qb,...
            TU,TUh,TQ,TQh,Tf,Tfh,TUhdir);
        
        
        % E_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/E.txt');
        % max(abs(col(E(:,:,1))-E_fort(:)))
        % stop
        
        % Cv_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/Cv.txt');
        % H_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/H.txt');
        %
        % max(abs(col(Cv(:,:,1)-TU(:,:,1))-Cv_fort(:)))
        % max(abs(col(H(:,:,1)-TUh(:,:,1))-H_fort(:)))
        % stop
        
        % UU_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/UU.txt');
        % UU0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/U0.txt');
        % LL_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/LL.txt');
        % LL0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/L0.txt');
        % max(abs(col(UU(:,:,1))-UU_fort(:)))
        % max(abs(col(UU0(:,1))-UU0_fort(:)))
        % max(abs(col(LL(:,:,1))-LL_fort(:)))
        % max(abs(col(LL0(:,1))-LL0_fort(:)))
        % stop
        % Impose boundary condition of the global problem
        [Df,Ef,Hf,Lf,Qf,fH,TUhf,TQhf,Tfhf] = hdg_bc_3D(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,T,infoFaces,refEl,boundaryNames,F,F_dir,u_tilde,gradU);
        
        %          Hf = zeros(size(Hf));
        
        %         %         % Checking
%                         residual = computeResiduals_local_analiticalGrad_3D...
%                             (F,F_dir,flipFace,u_tilde,u0,M,G,Cv,D,E,Edir,H,Hdir,force,...
%                             B,C,C_dir,L,P,Pb,Q,Qb,dt,invL,TU,TUh,TQ,TQh,Tf,Tfh,TUhdir,refEl)
%                 stop
        %                 resglob  = computeResiduals_global...
        %                     (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Qf,Hf,fH,F_dir,TUhf,TQhf,Tfhf);
        %                 sum(sum(resglob))
        %                 stop
        %         resglob_alt  = computeResiduals_global_alt...
        %             (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Hf,fH,...
        %             LL,LL0,UU,UU0,F_dir)
        %                stop
        % checkFortMapping
        % %
        % stop
        % Assembly
        [K,f] = hdg_Assembly_3D(refEl,F,F_dir,LL,LL0,UU,UU0,Df,Ef,Lf,Hf,Qf,fH, TUhf,TQhf,Tfhf);
       
%         figure, spy(K)
%         stop
        
        
%         
                debugFortran
                stop
        
        %         clear TU TUh TQ TQh Tf Tfh TUhdir TUhf TQhf Tfhf Cv H Hdir Hf
        %  matK= readCSRtxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/MatK.txt');
        % % %  matK = readCSRtxtParall('/home/giorgio/Dropbox/Fortran/MHDG/test/test_parall/MatK',2);
        % % %  f_fort = readCSRVectortxtParall('/home/giorgio/Dropbox/Fortran/MHDG/test/test_parall/rhs',2);
        % f_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/f.txt');
        % % % figure, spy(matK)
        % % % figure, spy(K)
        % disp(['Error in global matrix: ' num2str(max(max(abs(matK-K))))])
        % disp(['Error in global RHS: ' num2str((max(abs(f-f_fort))))])
        % stop
        
        %         Hf = zeros(size(Hf));
        %         [Kaux,faux] = hdg_Assembly(F,nv,U,U0,Df,Ef,Hf,fH,F_dir,nIntFaces);
        %         max(max(abs(K-Kaux)))
        %         stop
        %         condest(K)/1e9
        %         stop
        % face solution
        if savememory
            save('allvar.mat')
            clear
            load('allvar.mat','K','f','savememory')
        end
        sol = K\f;
        disp(['Cond(K) = ',num2str(condest(K),'%e')])
        
        if savememory
            save('allvar.mat','sol','-append')
            load('allvar.mat')
        end
%         clear K f
        
        
%         u_tilde = [sol(1:nvUnk); zeros(nvTot-nvUnk,1)];
        u_tilde = dumpNR*[sol(1:nvUnk); zeros(nvTot-nvUnk,1)]+(1-dumpNR)*u_tilde;
        
% u_tilde_fort =  readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/u_tilde.txt']);
% max(max(abs(u_tilde-u_tilde_fort)))
% stop
        % element by element solution
        u = computeElementByElementSolution_3D(ntor,u_tilde,UU,UU0,T,F,F_dir,refEl,refElTor );
        gradU = computeElementByElementGradient_3D(ntor,u_tilde,LL,LL0,T,F,F_dir,refEl,refElTor );
        %         clear LL LL0 UU UU0
        
%         
% u_fort =  readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/u.txt']);
% gradU_fort =  readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/q.txt']);
% max(max(abs(u-u_fort)))
% max(max(abs(gradU-gradU_fort)))
% stop        
        %% Check and fix problems in the solution
        if useThreshold~=0
            u = checksolution(u,X,T,refEl);
        end
        
        if drawsol
            up = cons2phys(transpose(reshape(u,[neq,numel(u)/neq])));
            figure(100+ifig)
            subplot(2,5,1),plotSolution(X,T,up(:,1),refEl,5);  title('Density')
            subplot(2,5,2), plotSolution(X,T,up(:,2),refEl,5); title('Parallel velocity')
            subplot(2,5,3),plotSolution(X,T,up(:,3),refEl,5);  title('Total energy for ions ')
            subplot(2,5,4), plotSolution(X,T,up(:,4),refEl,5); title('Total energy for electrons')
            subplot(2,5,5), plotSolution(X,T,up(:,5),refEl,5); title('Ions pressure HDG')
            subplot(2,5,6), plotSolution(X,T,up(:,6),refEl,5); title('Electrons pressure HDG')
            subplot(2,5,7), plotSolution(X,T,up(:,7),refEl,5); title('Ions temperature HDG')
            subplot(2,5,8), plotSolution(X,T,up(:,8),refEl,5); title('Electrons temperature HDG')
            subplot(2,5,9), plotSolution(X,T,up(:,9),refEl,5); title('Sound speed HDG')
            subplot(2,5,10), plotSolution(X,T,up(:,10),refEl,5); title('Mach')
            % % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
            drawnow
            
%             subplot(2,2,1),plotSolution(X,T,u(1:4:end-3),refEl); title('U1')
%             subplot(2,2,2),plotSolution(X,T,u(2:4:end-2),refEl); title('U2')
%             subplot(2,2,3),plotSolution(X,T,u(3:4:end-1),refEl); title('U3')
%             subplot(2,2,4),plotSolution(X,T,u(4:4:end),refEl);    title('U4')
%             drawnow
            
            
        end
        
        
        % Checking
        %         checkMapping(LL,LL0,UU,UU0,u_tilde,F,T,X);
        %         residual = computeResiduals_local...
        %             (Ne,Nv,nv,F,flipFace,u_tilde,u,M,G,Cv,D,E,Edir,H,Hdir,force,fdiff,...
        %             B,C,C_dir,L,P,Pb,Q,Qb,dt,u0,gradU);
        %         sum(sum(residual))
        %         stop
        
        
        
        %% Check if we got imaginary parts
        if ~isreal(u),error('Complex part detected'),end
        
        % check convergence for Newton Raphson
        errorNR = norm(u-u_iter)/sqrt(numel(u))/dumpNR;
        if errorNR<tolNR
            disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
            break
        elseif errorNR > 1e6
            error('Problem in the N-R procedure')
        else
            
            %             if iter>=7
            %                 up1 = cons2phys(transpose(reshape(u,[neq,numel(u)/neq])));
            %                 up0 = cons2phys(transpose(reshape(u_iter,[neq,numel(u_iter)/neq])));
            %                 udiff = up1-up0;
            %                 figure, plotSolution(X,T,udiff(:,2),refEl); title('Parallel velocity HDG')
            %             end
            
            u_iter = u;
            disp(['Error in Newton Raphson iteration '  num2str(iter) ': ' num2str(errorNR)])
            
            %         up = cons2phys(transpose(reshape(u,[neq,numel(u)/neq])));
            %         figure(100)
            %         subplot(2,4,1),plotSolution(X,T,up(:,1),refEl,5);  title('Density HDG')
            %         subplot(2,4,2), plotSolution(X,T,up(:,2),refEl,5); title('Parallel velocity HDG')
            %         subplot(2,4,3), plotSolution(X,T,up(:,3),refEl,5); title('Total energy HDG')
            %         subplot(2,4,4), plotSolution(X,T,up(:,4),refEl,5); title('Pressure HDG')
            %         subplot(2,4,5), plotSolution(X,T,up(:,5),refEl,5); title('Temperature HDG')
            %         subplot(2,4,6), plotSolution(X,T,up(:,6),refEl,5); title('Sound speed HDG')
            %         subplot(2,4,7), plotSolution(X,T,up(:,7),refEl,5); title('Mach HDG')
            % %         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
            %         drawnow
            %         stop
            %
            % figure(101)
            % subplot(1,3,1),plotSolution(X,T,u(1:3:end-2),refEl,5);  title('U1')
            % subplot(1,3,2), plotSolution(X,T,u(2:3:end-1),refEl,5); title('U2')
            % subplot(1,3,3), plotSolution(X,T,u(3:3:end),refEl,5);    title('U3')
            %         drawnow
        end
        
        
        
        % figure(10)
        % clf
        % subplot(1,3,1),plotSolution(X,T,u(1:3:end-2),refEl); title('Density HDG')
        % subplot(1,3,2),plotSolution(X,T,u(2:3:end-1)./u(1:3:end-2)./up(:,6),refEl); title('Parallel Mach HDG')
        % subplot(1,3,3),plotSolution(X,T,u(3:3:end)./u(1:3:end-2)./up(:,6),refEl); title('Energy HDG')
        % drawnow
        
    end
    
    disp(['Newton Raphson iterations: ' num2str(iter)])
    %     if iter == maxIter && iter~=1
    %         error('Problem: Newton Raphson does not converge')
    %     end
    
    % plot
    if ~mod(iStep,nPlot)
        converged_case = 0;
        saveSolutionWithParametersName;
        %        up = cons2phys(transpose(reshape(u,[neq,numel(u)/neq])));
        %         figure(100)
        %         subplot(2,4,1),plotSolution(X,T,up(:,1),refEl,5);  title('Density HDG')
        %         subplot(2,4,2), plotSolution(X,T,up(:,2),refEl,5); title('Parallel velocity HDG')
        %         subplot(2,4,3), plotSolution(X,T,up(:,3),refEl,5); title('Total energy HDG')
        %         subplot(2,4,4), plotSolution(X,T,up(:,4),refEl,5); title('Pressure HDG')
        %         subplot(2,4,5), plotSolution(X,T,up(:,5),refEl,5); title('Temperature HDG')
        %         subplot(2,4,6), plotSolution(X,T,up(:,6),refEl,5); title('Sound speed HDG')
        %         subplot(2,4,7), plotSolution(X,T,up(:,7),refEl,5); title('Mach HDG')
        %         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        %         drawnow
        
        
        % figure(101)
        % subplot(1,3,1),plotSolution(X,T,u(1:3:end-2),refEl,5);  title('U1')
        % subplot(1,3,2), plotSolution(X,T,u(2:3:end-1),refEl,5); title('U2')
        % subplot(1,3,3), plotSolution(X,T,u(3:3:end),refEl,5);    title('U3')
        %         drawnow
        
    end
    
    
    %     figure, plotSolution(X,T,u(1:2:end-1),refEl);hold on, plotMesh(X,T)
    %     figure, plotSolution(X,T,u0(1:2:end-1),refEl);hold on, plotMesh(X,T)
    %     stop
    
    
    % check convergence for the time iterations
    errorTime = norm(u-u0)/sqrt(numel(u));
    if errorTime < tolTime
        if strcmpi(problem,'NS_time')
%             disp(['Error in time ' num2str(errorTime)])
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
%         disp(['Error in time step '  num2str(iStep) ': ' num2str(errorTime)])
%         disp(['Time step: ' num2str(dt)])
%         disp(['Simulation time: ' num2str(time)])
%         disp(' ')
    end
    
    
end

toc
ifig = ifig+1;
% store_simulation_parameters;
% converged_case = 1;
% saveSolutionWithParametersName;

%  disp(['Cond(K) = ',num2str(condest(K),'%e')])


up = cons2phys(transpose(reshape(u,[neq,numel(u)/neq])));

% Plot solution
ur = transpose(reshape(u,[neq,numel(u)/neq]));
disp('Plotting...')
for itor = 1:ntpos
    %     upol = extractSolutionInAPoloidalPlane(u,T,refElPol,refElTor,itor);
    upol = extractSolutionInAtGivenTheta(ur(:,1),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - U1'])
    upol = extractSolutionInAtGivenTheta(ur(:,2),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - U2'])    
    upol = extractSolutionInAtGivenTheta(ur(:,3),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - U3'])    
    upol = extractSolutionInAtGivenTheta(ur(:,4),T,refEl,refElTor,tpos(itor));
    figure
    plotSolution(X,T,upol,refEl)
    title(['Plane ' num2str(itor) ' - U4'])    
    
end



% disp('Plotting analytical sol...')
% uan = initializeSolutionToAnalyticalSol_3D(X,T,refEl);
% uan = transpose(reshape(uan,[neq,numel(uan)/neq]));
% 
% for itor = 1:ntpos
%     %     upol = extractSolutionInAPoloidalPlane(u,T,refElPol,refElTor,itor);
%     upol = extractSolutionInAtGivenTheta(uan(:,1),T,refEl,refElTor,tpos(itor));
%     figure
%     plotSolution(X,T,upol,refEl)
%     title(['Plane ' num2str(itor) ' - U1 - anal'])
%     upol = extractSolutionInAtGivenTheta(uan(:,2),T,refEl,refElTor,tpos(itor));
%     figure
%     plotSolution(X,T,upol,refEl)
%     title(['Plane ' num2str(itor) ' - U2 - anal'])    
%     upol = extractSolutionInAtGivenTheta(uan(:,3),T,refEl,refElTor,tpos(itor));
%     figure
%     plotSolution(X,T,upol,refEl)
%     title(['Plane ' num2str(itor) ' - U3 - anal'])    
%     upol = extractSolutionInAtGivenTheta(uan(:,4),T,refEl,refElTor,tpos(itor));
%     figure
%     plotSolution(X,T,upol,refEl)
%     title(['Plane ' num2str(itor) ' - U4 - anal'])    
%     
% end

%% L2 error
disp('Computing error...')
err = computeError_3D(u,X,T,theta,refEl,refElTor,ntor);
disp(['Error U1 = ', num2str(err(1),10)]);
disp(['Error U2 = ', num2str(err(2),10)]);
disp(['Error U3 = ', num2str(err(3),10)]);
disp(['Error U4 = ', num2str(err(4),10)]);

disp(' ')


%% Analytical solution
% uan = analyticalSolution(X);


% uan = analyticalSolution(X);

% figure, plotSolution(X,T,up(:,1),refEl); title('Density HDG')
% figure, plotSolution(X,T,up(:,2),refEl); title('Parallel velocity HDG')
% figure, plotSolution(X,T,up(:,3),refEl); title('Total energy HDG')
% figure, plotSolution(X,T,up(:,4),refEl); title('Pressure HDG')
% figure, plotSolution(X,T,up(:,5),refEl); title('Temperature HDG')
% figure, plotSolution(X,T,up(:,6),refEl); title('Sound speed HDG')
% figure, plotSolution(X,T,up(:,7),refEl); title('Mach HDG')


% figure(100+ifig)
% subplot(2,5,1),plotSolution(X,T,up(:,1),refEl,5);  title('Density')
% subplot(2,5,2), plotSolution(X,T,up(:,2),refEl,5); title('Parallel velocity')
% subplot(2,5,3),plotSolution(X,T,up(:,3),refEl,5);  title('Total energy for ions ')
% subplot(2,5,4), plotSolution(X,T,up(:,4),refEl,5); title('Total energy for electrons')
% subplot(2,5,5), plotSolution(X,T,up(:,5),refEl,5); title('Ions pressure HDG')
% subplot(2,5,6), plotSolution(X,T,up(:,6),refEl,5); title('Electrons pressure HDG')
% subplot(2,5,7), plotSolution(X,T,up(:,7),refEl,5); title('Ions temperature HDG')
% subplot(2,5,8), plotSolution(X,T,up(:,8),refEl,5); title('Electrons temperature HDG')
% subplot(2,5,9), plotSolution(X,T,up(:,9),refEl,5); title('Sound speed HDG')
% subplot(2,5,10), plotSolution(X,T,up(:,10),refEl,5); title('Mach')
% % % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% drawnow





% figure
% subplot(2,2,1),plotSolution(X,T,uan(:,1),refEl); title('U1 analytical')
% subplot(2,2,2),plotSolution(X,T,uan(:,2),refEl); title('U2 analytical')
% subplot(2,2,3),plotSolution(X,T,uan(:,3),refEl); title('U3 analytical')
% subplot(2,2,4),plotSolution(X,T,uan(:,4),refEl); title('U4 analytical')



% figure
% subplot(2,2,1),plotSolution(X,T,u(1:4:end-3),refEl,5); title('U1')
% subplot(2,2,2),plotSolution(X,T,u(2:4:end-2),refEl,5); title('U2')
% subplot(2,2,3),plotSolution(X,T,u(3:4:end-1),refEl,5); title('U3')
% subplot(2,2,4),plotSolution(X,T,u(4:4:end),refEl,5);    title('U4')



% L2err = computeErrorAnalyticSol(X,T,u,refEl);
% disp(['Error U1: ' num2str(L2err(1))])
% disp(['Error U2: ' num2str(L2err(2))])
% disp(['Error U3: ' num2str(L2err(3))])
% disp(['Error U4: ' num2str(L2err(4))])
% subplot(2,2,4),plotSolution(X,T,uan(:,2)./uan(:,1),refEl); title('Parallel velocity analytical')



% umat = u;
% HDF5load('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/Sol_Circ_3_P3_Diff.20000E+01.h5')
% u = u';


% if driftvel
% %     plotDiamagDrift(X,T)
%     if (testcase.n >= 50 && testcase.n<60)
%         diamagvel = sqrt(Magnetic.dvxnodes.^2+Magnetic.dvynodes.^2);
%         diamagvel = col(diamagvel);
%     else
%         [~,~,dvx,dvy] = defineMagneticField(X);
%         diamagvel = sqrt(dvx.^2+dvy.^2);
%     end
%     plotSolution(X,T,diamagvel,refEl); title('Diamagnetic velocity')
%
% end
%% plot solution on a line
% figure
% line = [linspace(2.5,3,100)'/Mesh.lscale,zeros(100,1)];
% uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
% semilogy(line(:,1)*Mesh.lscale,uplot),ylim([1e-2 1e0])
% hold on
% lambda = sqrt(43*2/1.39e5);
% plot(line(:,1)*Mesh.lscale,line(1,1)*Mesh.lscale*lambda+1-line(:,1)*Mesh.lscale*lambda,'r')
%
% L2err = computeErrorAnalyticSol(X,T,u,refEl)
% plotError(X,T,u,refEl)

% u0 = u; save u_WEST_h0.02_Diff0.04_P4.mat u0


% plot
% plotVelocity(X,T,u)
% plotPressure(Xp,Tp,p)
% figure
% subplot(2,2,1),plotDiscontinuousSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
% subplot(2,2,2),plotDiscontinuousSolution(X,T,u(2:2:end)./u(1:2:end-1),refEl); title('Parallel velocity HDG')
% subplot(2,2,3),plotSolution(X,T,uan(:,1),refEl); title('Density analytical')
% subplot(2,2,4),plotSolution(X,T,uan(:,2)./uan(:,1),refEl); title('Parallel velocity analytical')


% figure
% subplot(2,2,1),plotDiscontinuousSolution(X,T,u(1:2:end-1),refEl); title('U(1) HDG')
% subplot(2,2,2),plotDiscontinuousSolution(X,T,u(2:2:end),refEl); title('U(2) HDG')
% subplot(2,2,3),plotSolution(X,T,uan(:,1),refEl); title('U(1) analytical')
% subplot(2,2,4),plotSolution(X,T,uan(:,2),refEl); title('U(2) velocity analytical')
%
% L2err = computeErrorAnalyticSol(X,T,u,refEl);
% disp(['L2 error density: ', num2str(L2err(1))])
% disp(['L2 error parallel velocity: ', num2str(L2err(2))])

% figure, plotDiscontinuousSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
% figure, plotDiscontinuousSolution(X,T,u(2:2:end)./u(1:2:end-1),refEl); title('Parallel velocity HDG')


% figure, plotSolution(X,T,uan(:,1),refEl); title('Density analytical')
% figure, plotSolution(X,T,uan(:,2)./uan(:,1),refEl); title('Parallel velocity analytical')

% % streamlines
% u_cont = createContinuousSolution_vec(X,T,u);
% [sx,sy] = defineStartingPoendints;
% plotStreamLines(X,T,u_cont,100,sx,sy)

% % postprocess
% L_gradv = calculateElementByElementVelocityGradient(u_tilde,ro,F,L,L0,Lro,Lf,U);
% [u_star u_int] = hdg_PostprocessSolution(L_gradv,u,Kp,Btp,int_Np,refElv_star,refEl,Ne);
% %
% % streamlines for the postprocessed velocity
% u_star_cont = createContinuousSolution_vec(X_star,T_star,u_star);
% [sx,sy] = defineStartingPoints;
% plotStreamLines(X_star,T_star,u_star_cont,100,sx,sy)
%
% % error
% error = computeErrorVelocity(M_star,u_star,u_int);
% figure, plotErrorMap(X,T,error.^0.5)
%

% % exact error
% [errorv_ex egv_ex errorp_ex egp_ex] = computeErrorAnalyticSol(Re,u,X,T,refEl,p,refElp);
% figure, plotErrorMap(X,T,errorv_ex) ,title('Error analytical')

% refSol = load(['Cavity_Re' num2str(Re) '.mat']);
% [error eg domArea] = calculateL2ErrorTwoSolutionNested(X,T,u,refEl,refSol.X,refSol.T,refSol.u,refSol.refEl);
% figure, plotErrorMap(X,T,error),title('Error nested')
% [error_post eg_post domArea] = calculateL2ErrorTwoSolutionNested(X_star,T_star,u_star,refElv_star,...
%     refSol.X,refSol.T,refSol.u,refSol.refEl);
% figure, plotErrorMap(X,T,error_post),title('Error nested post')


% disp(['Error velocity exact solution: ' num2str(egv_ex)])
% disp(['Error pressure exact solution: ' num2str(egp_ex)])
%
% plotStreamlinesWithStreamFunction(X,T,u_cont,refEl,boundary,100)
% %
% plotStreamlinesWithStreamFunction(X_star,T_star,u_star_cont,refElv_star,boundary_star,100)
%
% save data
% saveData

% % analytical solution
% [u_an_cont p_an] = analyticalSolution(X);
% u_an = zeros(numel(u_an_cont),1);
% u_an(1:2:end-1) = u_an_cont(:,1);
% u_an(2:2:end) = u_an_cont(:,2);
% plotStreamlinesWithStreamFunction(X,T,u_an,refEl,boundary,100)
% figure, plotMesh(X,T); plotSolution(X,T,p_an,refEl); axis equal

