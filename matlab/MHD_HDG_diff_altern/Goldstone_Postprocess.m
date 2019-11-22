% Postprocess routine

clear
% close all

global testcase Mesh

lscale = 1.901*1e-3;
kT = 1;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
linespec = {'b-','b--','r-','r--','k-','k--','m-','m--'};

figure

points = 1:9;
D = 0.1;
Diffvec = D*1./sqrt(10).^(points-1);

%% ************************************************************************
%
%%                          Cases without Drift 
%
%% ************************************************************************
disp('  ')
disp('Cases without drift')
disp('  ')

%% Diffusion 0.1 - 0.00316

% Load mesh
meshName = 'Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

% Initialize
lambda_n = zeros(size(points,2),1);
lambda_n_est = lambda_n;
leg = cell(numel(points),1);
for idiff = 1:4
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_'...
        'refSep0.01_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    np = 100;
    line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.85/lscale) == min(abs(line(:,1)-2.85/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.86/lscale) == min(abs(line(:,1)-2.86/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end


%% Diffusion 0.001 - 0.000316

% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

for idiff = 5:6
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.'...
        'mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    np = 100;
    line = [linspace(2.5,2.8,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.765/lscale) == min(abs(line(:,1)-2.765/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.775/lscale) == min(abs(line(:,1)-2.775/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end

%% Diffusion 0.0001

% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P6.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

for idiff = 7:7
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_'...
        'refSep0.008_P6.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1_LocDiffPts2.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    np = 100;
    line = [linspace(2.5,2.77,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.755/lscale) == min(abs(line(:,1)-2.755/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.76/lscale) == min(abs(line(:,1)-2.76/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end


% legend(leg{1:7})


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
meshName = 'Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

% Initialize
lambda_n_drift = zeros(size(points,2),1);
lambda_n_est_drift = lambda_n_drift;
leg_drift = cell(numel(points),1);
figure
for idiff = 1:3
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_'...
        'refSep0.01_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1_Drift_LocDiffPts2_UseTh1e-06.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    np = 100;
    line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.85/lscale) == min(abs(line(:,1)-2.85/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.86/lscale) == min(abs(line(:,1)-2.86/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est_drift(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est_drift(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end



%% Diffusion 0.00316 - 0.001

% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';

for idiff = 4:5
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_'...
        'refSep0.008_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1_Drift_LocDiffPts2_UseTh1e-06.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %

    np = 100;
    line = [linspace(2.5,2.8,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.765/lscale) == min(abs(line(:,1)-2.765/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.775/lscale) == min(abs(line(:,1)-2.775/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est_drift(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est_drift(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end

%% Diffusion 0.00031623 - 0.0001

% Load mesh
meshName = 'Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';
for idiff = 6:7
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = './Saves/';
    solName = ['u_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_'...
        'DoubleSol_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1_Drift_SC2_LocDiffPts2.mat'];     
    load([path2Mesh meshName])
    load([path2Sol solName])
    
    % Apply translation for axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - 0.5*(min(X(:,1))+max(X(:,1)))+2;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
        
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %

    np = 100;
    line = [linspace(2.5,2.77,np)'/lscale,zeros(np,1)];    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.755/lscale) == min(abs(line(:,1)-2.755/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.76/lscale) == min(abs(line(:,1)-2.76/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est_drift(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est_drift(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
end


figure(10),loglog(Diffvec(lambda_n~=0),lambda_n(lambda_n~=0)/lscale,'k-o');
hold on,   loglog(Diffvec(lambda_n_drift~=0),lambda_n_drift(lambda_n_drift~=0)/lscale,'k-*');
hold on,   loglog(Diffvec(lambda_n~=0),31.6227766016838*sqrt(Diffvec(lambda_n~=0)),'k--')
legend('Without drift', 'With drift','Slope 0.5 reference','location','southeast')
xlabel('Diffusion')
ylabel('$\lambda_n (\rho_{Larm.})$','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'GoldstoneBW')
