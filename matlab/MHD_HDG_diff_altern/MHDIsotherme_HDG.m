% HDG method for the Isothermal Euler equations
clear
clear global
close all
home

global kT testcase Mesh refEl neq diff_n diff_u Magnetic coreflux axisym
global driftvel  shockcapt localDiffPoints limitRhoMin useThreshold dirichlet_weak
kT = 25;%12;
% kT =  0.024;
neq = 2;   % n of equations (2 for MHD isothermal 2D)
axisym = 0;
driftvel = 0;
shockcapt = 0; % 0- no sh 1-elementwise 2-linear interp
localDiffPoints = 0;
limitRhoMin = 0;
useThreshold = 0;
MultAddDiff = 1;

%% Define problem
% problem = 'NS_time';
problem = 'NS_steady';
% problem = 'NS_pseudotime';
% problem = 'S_steady';

%% Test case
testcase.n = 9;
testcase.wavex = 1;
testcase.wavey = 1;
testcase.xc = -0.5;
testcase.yc = -0.5;
setneum = 0; % Neumann BC
diff_n = 0;%3.8;
diff_u = 0;%3.8;
coreflux = [1;0];
dirichlet_weak = 0;
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

%% Initialization
if restart
    %     load('../MHD_HDG_diff_altern/Saves/u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_0.0031623_tau1_UseTh1e-06.mat','u0','u_tilde');
    load('../MHD_HDG_diff_altern/Saves/u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat_Diffusion_0.1_tau1.mat','u0','u_tilde');
    
%         u0 = projectFieldDifferentP(u0,4,refEl.degree,size(T,1));
else
    % u0 = initializeSolution(X,T,F_dir,refEl);
    u0 = initializeSolutionToAnalyticalSol(X,T);
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

if shockcapt
    shock_st = shockCaptPrep(refEl);
end


%% Set dimensions
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
h = 1;

%% Precalculated matrices
disp('Precalculated matrices...')
[M,B,C,C_dir,G,D,E,Df,Ef,Edir,invL,force,Lf,...
    L,P,Pb,Q,Qb,Qf] =...
    hdg_PrecalculatedMatrices(X,T,flipFace,refEl,tau,F,F_dir);


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
Cv = zeros(neq*Nv,neq*Nv,Ne);
H = zeros(neq*Nv,nf*neq*nv,Ne);
Hf = zeros(nf*neq*nv,nf*neq*nv,Ne);
Hdir = zeros(neq*Nv,Ne);
Hfdir = zeros(nf*neq*nv,Ne);
fH = Hfdir;


%% check the problem type
checkProblem

%% Extract u tilde
if ~exist('u_tilde','var')
    u_tilde = extractFaceSolution(u0,F,F_dir,refEl);
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
    
    ifig=0;
    %% Newton Raphson iterations
    for iter = 1:maxIter
        
        ifig = ifig+1;
        
        % Calculate convection matrices (only if NS)
        if strcmp(problem,'NS_time') || strcmp(problem,'NS_steady') || strcmp(problem,'NS_pseudotime')
            [Cv,H,Hdir,Hf] = ...
                hdg_ConvectionMatrices(X,T,F,flipFace,refEl,u_iter,u_tilde,F_dir);
        end

%                 residual = computeResiduals_local_analiticalGrad...
%                     (Ne,Nv,nv,F,flipFace,u_tilde,u0,M,G,Cv,D,E,Edir,H,Hdir,force,...
%                     B,C,C_dir,L,P,Pb,Q,Qb,dt,invL)
% stop
        % Get the mapping
        [LL,LL0,UU,UU0] = hdg_Mapping(flipFace,Nv,nv,M,G,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
            B,C,C_dir,invL,P+P_locDiffPts+P_sc+P_limRho,Pb,Q+Q_locDiffPts+Q_sc+Q_limRho,Qb);
        
        
        
%         checkMapping(LL,LL0,UU,UU0,u_tilde,F,T,X)
        % Impose boundary condition of the global problem
        [Df,Ef,Hf,Lf,Qf,fH] = hdg_bc(Df,Ef,Hf,Lf,Qf,X,T,infoFaces,refEl,boundaryNames,F,u_tilde);
        
        %          Hf = zeros(size(Hf));
        
        %         % Checking

%         stop
        %         resglob  = computeResiduals_global...
        %             (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Hf,fH,F_dir);
        %         sum(sum(resglob))
        %         resglob_alt  = computeResiduals_global_alt...
        %             (Ne,Nf,Nv,nv,F,u_tilde,u_iter,Df,Ef,Lf,Hf,fH,...
        %             LL,LL0,UU,UU0,F_dir)
        %                stop
% checkFortMapping 
% % 
% stop
        % Assembly
        [K,f] = hdg_Assembly(F,nv,LL,LL0,UU,UU0,Df,Ef,Lf+Lf_locDiffPts+Lf_sc+Lf_limRho,...
            Hf,Qf,fH,F_dir,nIntFaces);

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
        u = calculateElementByElementSolution(u_tilde,F,UU,UU0);
        gradU = calculateElementByElementGradient(u_tilde,F,LL,LL0);
 
% figure(ifig)        
% subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
% subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
% drawnow        
        % Checking
        %         checkMapping(LL,LL0,UU,UU0,u_tilde,F,T,X);
        %         residual = computeResiduals_local...
        %             (Ne,Nv,nv,F,flipFace,u_tilde,u,M,G,Cv,D,E,Edir,H,Hdir,force,fdiff,...
        %             B,C,C_dir,L,P,Pb,Q,Qb,dt,u0,gradU);
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


%% Analytical solution
% uan = reshape(analyticalSolution(X),size(X));
% subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
% subplot(1,2,2),plotSolution(X,T,u(2:2:end)./(u(1:2:end-1)+eps)/sqrt(kT),refEl); title('Parallel Mach HDG')
% plotParallelVelocityVector(X,T,u(2:2:end)./(u(1:2:end-1)+eps),10)


% uan = reshape(analyticalSolution(X),size(X));

figure
subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl,5); title('Density HDG')
subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl,5); title('Parallel velocity HDG')
% subplot(2,2,3),plotSolution(X,T,uan(:,1),refEl); title('Density analytical')
% subplot(2,2,4),plotSolution(X,T,uan(:,2)./uan(:,1),refEl); title('Parallel velocity analytical')




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
L2err = computeErrorAnalyticSol(X,T,u,refEl)
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

