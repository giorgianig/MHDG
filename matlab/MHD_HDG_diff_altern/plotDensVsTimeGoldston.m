% plot density vs time Gold
clear
close all

path2Mesh = '../../Meshes/Meshes_2D/';


figure
% Diff 3e-3
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-3/';  
solName = 'u_Circle_ONION_2_P6_Diff.31623E-02';
nsol = 76;
dt = 1000;
minrho = zeros(nsol,1);
time = dt*(1:nsol);
for isol = 1:nsol
    
    load([path2Sol,solName,'_',num2str(isol,'%04d'),'.mat'])
    minrho(isol) = min(u(1:2:end-1));
end
semilogy(time,minrho)
hold on


% Diff 1e-3
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-3/';  
solName = 'u_Circle_ONION_2_P6_Diff.10000E-02';
nsol = 424;
dt = 100;
minrho = zeros(nsol,1);
time = dt*(1:nsol);
for isol = 1:nsol
    
    load([path2Sol,solName,'_',num2str(isol,'%04d'),'.mat'])
    minrho(isol) = min(u(1:2:end-1));
end
semilogy(time,minrho)
hold on



% Diff 3e-4
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D3e-4/';  
solName = 'u_Circle_ONION_2_P6_Diff.31623E-03';
nsol = 448;
dt = 100;
minrho = zeros(nsol,1);
time = dt*(1:nsol);
for isol = 1:nsol
    
    load([path2Sol,solName,'_',num2str(isol,'%04d'),'.mat'])
    minrho(isol) = min(u(1:2:end-1));
end
semilogy(time,minrho)
hold on



% Diff 1e-4
path2Sol = '/home/giorgio/Goldston/Matlab/YesDrift/D1e-4/';  
solName = 'u_Circle_ONION_2_P10_Diff.10000E-03';
nsol = 310;
dt = 100;
minrho = zeros(nsol,1);
time = dt*(1:nsol);
for isol = 1:nsol
    
    load([path2Sol,solName,'_',num2str(isol,'%04d'),'.mat'])
    minrho(isol) = min(u(1:2:end-1));
end
semilogy(time,minrho)
hold on




legend('3e-3','1e-3','3e-4','1e-4')
grid on