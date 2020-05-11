% Postprocess routine

clear
 close all

global testcase Mesh

% lscale = 1.901*1e-3;
lscale = 1.*1e-3;
kT = 25;
testcase.n = 60;
axisym = 1;
fref = 7.28e6; % reference frequency
linespec = {'b-','b--','r-','r--','k-','k--','m-','m--'};

ndiv = 6;
np = 1000;
points = 1:9;





theta = linspace(0,2*pi,ndiv);
D = 0.1;
Diffvec = D*1./sqrt(10).^(points-1);


path2Mesh = '../../Meshes/Meshes_2D/';
% Initialize
lambda_n = zeros(size(points,2),1);
lambda_n_est = lambda_n;
leg = cell(numel(points),1);


%% ************************************************************************
%
%%                          Cases without Drift
%
%% ************************************************************************
disp('  ')
disp('Cases without drift')
disp('  ')

figure

% 

%% Diffusion 0.1 

% Load mesh
meshName = 'Circle_ONION_1_P6.mat';

for idiff = 1:1
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_1_P6_Diff.10000E+00.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.34/lscale) == min(abs(line(:,1)-0.34/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.35/lscale) == min(abs(line(:,1)-0.35/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% Diffusion 0.031 
% Load mesh
meshName = 'Circle_ONION_1_P6.mat';
for idiff = 2:2
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_1_P6_Diff.31623E-01.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));

    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.34/lscale) == min(abs(line(:,1)-0.34/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.35/lscale) == min(abs(line(:,1)-0.35/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% Diffusion 0.01 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff = 3:3
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_2_P6_Diff.10000E-01.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end




%% Diffusion 0.0031 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff =4:4
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_2_P6_Diff.31623E-02.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% Diffusion 0.001 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff = 5:5
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_2_P6_Diff.10000E-02.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end







%% Diffusion 0.00031 
% Load mesh
meshName = 'Circle_ONION_3_P6.mat';
for idiff = 6:6
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_3_P6_Diff.31623E-03.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.293/lscale) == min(abs(line(:,1)-0.293/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.294/lscale) == min(abs(line(:,1)-0.294/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
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
meshName = 'Circle_ONION_3_P6.mat';
for idiff = 7:7
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/NoDrift/';  
    solName = 'u_Circle_ONION_3_P6_Diff.10000E-03.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.292/lscale) == min(abs(line(:,1)-0.292/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.293/lscale) == min(abs(line(:,1)-0.293/lscale)));
    ien = ien(1);
    lambda_n(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% ************************************************************************

%%                          Cases with Drift

%% ************************************************************************

disp('  ')
disp('Cases with drift')
disp('  ')

figure

%% Diffusion 0.1 

% Load mesh
meshName = 'Circle_ONION_1_P6.mat';

for idiff = 1:1
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';  
    solName = 'u_Circle_ONION_1_P6_Diff.10000E+00.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.34/lscale) == min(abs(line(:,1)-0.34/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.35/lscale) == min(abs(line(:,1)-0.35/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% Diffusion 0.031 
% Load mesh
meshName = 'Circle_ONION_1_P6.mat';
for idiff = 2:2
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';  
    solName = 'u_Circle_ONION_1_P6_Diff.31623E-01.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));

    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.34/lscale) == min(abs(line(:,1)-0.34/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.35/lscale) == min(abs(line(:,1)-0.35/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



%% Diffusion 0.01 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff = 3:3
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/';  
    solName = 'u_Circle_ONION_2_P6_Diff.10000E-01.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end




%% Diffusion 0.0031 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff =4:4
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-3/';  
    solName = 'u_Circle_ONION_2_P6_Diff.31623E-02_0076.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end



% Diffusion 0.001 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff = 5:5
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-3/';  
    solName = 'u_Circle_ONION_2_P6_Diff.10000E-02_0424.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end





% 
% 
%% Diffusion 0.00031 
% Load mesh
meshName = 'Circle_ONION_2_P6.mat';
for idiff = 6:6
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-4/';  
    solName = 'u_Circle_ONION_2_P6_Diff.31623E-03_0448.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.32/lscale) == min(abs(line(:,1)-0.32/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.33/lscale) == min(abs(line(:,1)-0.33/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end





%% Diffusion 0.0001 
% Load mesh
meshName = 'Circle_ONION_2_P10.mat';
for idiff = 7:7
    
    disp(['Diffusion ', num2str(Diffvec(idiff))])
    leg{idiff} = ['Diffusion ', num2str(Diffvec(idiff))];
    
    %% Load solution
    path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-4/';  
    solName = 'u_Circle_ONION_2_P10_Diff.10000E-03_0310.mat';
    load([path2Mesh meshName])
    load([path2Sol solName])
    maxx = max(X(:,1));
    X = X/lscale;
    
    Mesh.X = X;
    Mesh.T = T;
    Mesh.lscale = lscale;
    
    
    refEl = createReferenceElement(1,size(T,2),[]);
    %
    
    line = [linspace(0.2,maxx,np)'/lscale,zeros(np,1)];
    
    
    line_rot = zeros(size(line,1)*ndiv,size(line,2));
    % rotation
    for it = 1:ndiv
        A = [cos(theta(it)) sin(theta(it));-sin(theta(it)) cos(theta(it))];
        line_rot((it-1)*np+(1:np),:) = transpose(A*line');
    end
    
    %     figure
    %     plotMesh(X,T)
    %     hold on
    %     plot(line_rot(:,1),line_rot(:,2),'r*')
    %     stop
    uplot = evalDGapproximationAtPoints(line_rot,u(1:2:end-1),X,T,refEl);
    uplot = reshape(uplot,np,ndiv);
    uplot = sum(uplot,2)/ndiv;

    semilogy(line(:,1)*lscale,uplot,'r-'),grid on
    hold on
    line(uplot==0,:) = [];
    uplot(uplot==0,:) = [];
    
    ist = find(abs(line(:,1)-0.292/lscale) == min(abs(line(:,1)-0.292/lscale)));
    ist = ist(1);
    ien = find(abs(line(:,1)-0.293/lscale) == min(abs(line(:,1)-0.293/lscale)));
    ien = ien(1);
    lambda_n_drift(idiff) = -(line(ien,1)-line(ist,1))*lscale/log10(uplot(ien)/uplot(ist));
    Diff = Diffvec(idiff)*fref*(lscale)^2;
    vref = lscale*fref;
    lambda_n_est(idiff) = sqrt(3/4*43*Diff/vref);
    disp(['Lambda_n from graphic: ' num2str(lambda_n_drift(idiff))])
    disp(['Lambda_n estimated: ' num2str(lambda_n_est(idiff))])
    disp(['Minimum density: ', num2str(min(uplot))])
    disp(' ')
    
end


figure(10),loglog(Diffvec(lambda_n~=0),lambda_n(lambda_n~=0)/lscale,'ko','Markersize',12);
hold on,   loglog(Diffvec(lambda_n_drift~=0),lambda_n_drift(lambda_n_drift~=0)/lscale,'k*','Markersize',12);
hold on,   loglog(Diffvec(lambda_n~=0),111.6227766016838*sqrt(Diffvec(lambda_n~=0)),'k--')
legend({'Without drift', 'With drift','$\propto \sqrt{D}$'},'interpreter','latex','location','southeast')
xlabel('$D$','interpreter','latex')
ylabel('$\lambda_n (\rho_{L})$','interpreter','latex')
grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'GoldstoneBW')
