% Postprocess routine

clear
% close all

global testcase Mesh





D = 0.1;
lscale = 1.901*1e-3;
kT = 1;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
points = 1:6;
% points = 1:2;
% Diffvec = D*1./sqrt(10).^(points-1);
Diffvec = 0.01;
%% Load mesh
% meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat';
% meshName = 'Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat';
% meshName = 'Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat';
meshName = 'Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat';
path2Mesh = '../../Meshes/Meshes_2D/';
lambda_n = zeros(size(points,1),1);
lambda_n_est = lambda_n;
linespec = {'b-','b--','r-','r--','k-','k--','m-','m--'};
leg = cell(numel(points),1);
% figure
for idiff = 1:numel(Diffvec)
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    % meshName = 'West_h0.04_refCorn_P8_fixed.mat';
    % meshName = 'West_SmallHole_h0.04_refCorn_P8_fixed.mat';
    % meshName = 'West_NoHole_h0.04_refCorn_P8_fixed.mat';
    % path2Mesh = '../../Meshes/Meshes_2D/WEST/';
    
    %% Load solution
    path2Sol = './Saves/';
    % solName = ['u_LIM_h01_refCorn_diff' num2str(D) '_P10.mat'];
    % solName = ['u_West_Case51_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
    % solName = ['u_West_SmallHole_Case52_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
    % solName = ['u_West_NoHole_Case53_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
%     solName = ['u_Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1.mat'];
%     solName = ['u_Circle_LIM_InfThin_h0.02EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P10.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1.mat'];
%     solName = 'u_Circle_LIM_InfThin_h0.05EXT_0.2INT_refCorn0.001_refBound0.02_refSep0.01_P5.mat_Diffusion_0.0031623_tau1LimitMinRho1e-6_NONCONV.mat';
solName = 'u_Circle_LIM_InfThin_SmallSOL_refCorn0.001_refSol0.005_refSep0.008_P10.mat_Diffusion_0.01_tau1.mat';
    %% Physical parameters
    
    
    
    %% Working...
    
    load([path2Mesh meshName])
    load([path2Sol solName])
    %% Apply translation if axisymmetric case
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
    % subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
    % subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
    % figure,plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG'),axis off
    % figure,plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG'),axis off
    
    % plotParallelVelocityVector(X,T,u(2:2:end)./u(1:2:end-1),10)
    
    % hold on
    % figure(2)
    hold on
    np = 1000;
    line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
    % line = [linspace(1.8,3,100)'/lscale,zeros(100,1)];
    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
%     hold on
    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.77/lscale) == min(abs(line(:,1)-2.77/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-2.78/lscale) == min(abs(line(:,1)-2.78/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))]) 
    disp(['Minimum density: ', num2str(min(uplot))])
end
legend(leg)

figure(10),plot(Diffvec,lambda_n/lscale,'b-o');

points = 1:5;
% points = 1:2;
lambda_n_df = zeros(size(points,1),1);
lambda_n_est = lambda_n_df;
Diffvec = D*1./sqrt(10).^(points-1);
leg = cell(numel(points),1);
figure
for idiff = 1:numel(Diffvec)
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    % meshName = 'West_h0.04_refCorn_P8_fixed.mat';
    % meshName = 'West_SmallHole_h0.04_refCorn_P8_fixed.mat';
    % meshName = 'West_NoHole_h0.04_refCorn_P8_fixed.mat';
    % path2Mesh = '../../Meshes/Meshes_2D/WEST/';
    
    %% Load solution
    path2Sol = './Saves/';
    % solName = ['u_LIM_h01_refCorn_diff' num2str(D) '_P10.mat'];
    % solName = ['u_West_Case51_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
    % solName = ['u_West_SmallHole_Case52_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
    % solName = ['u_West_NoHole_Case53_h0.04_refCorn_diff' num2str(D) '_P8_DiamagDrifg.mat'];
    solName = ['u_Circle_LIM_0.1_RefSepar0.01_RefCorn_0.001_RefOUT_0.02_P8.mat_Diffusion_' num2str(Diffvec(idiff)) '_tau1diffCorn1e-1_Drift.mat'];
    %% Physical parameters
    
    
    
    %% Working...
    
    load([path2Mesh meshName])
    load([path2Sol solName])
    %% Apply translation if axisymmetric case
    if axisym && min(X(:,1))<0
        % apply translation in x to get min(X(:,1)) = 1
        X(:,1) = X(:,1) - min(X(:,1))+1;
    end
    
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    u = u0; clear u0
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    % subplot(1,2,1),plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG')
    % subplot(1,2,2),plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG')
    % figure,plotSolution(X,T,u(1:2:end-1),refEl); title('Density HDG'),axis off
    % figure,plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl); title('Parallel Mach HDG'),axis off
    
    % plotParallelVelocityVector(X,T,u(2:2:end)./u(1:2:end-1),10)
    
    % hold on
    % figure(2)
    % hold on
    np = 1000;
    line = [linspace(2.5,3,np)'/lscale,zeros(np,1)];
    % line = [linspace(1.8,3,100)'/lscale,zeros(100,1)];
    
    uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
    semilogy(line(:,1)*lscale,uplot,linespec{idiff}),grid on
    hold on
    line(uplot==0) = [];
    uplot(uplot==0) = [];
    
    ist = find(abs(line(:,1)-2.75/lscale) == min(abs(line(:,1)-2.75/lscale)));
    ien = find(abs(line(:,1)-2.755/lscale) == min(abs(line(:,1)-2.755/lscale)));
    lambda_n_df(idiff) = -(line(ien,1)-line(ist,1))*lscale/log(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_df(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))]) 
end
legend(leg)

figure(10),hold on,plot(Diffvec,lambda_n_df/lscale,'r--o');
legend('Without drift', 'With drift')
xlabel('Diffusion')
ylabel('$\lambda_n$','interpreter','latex')