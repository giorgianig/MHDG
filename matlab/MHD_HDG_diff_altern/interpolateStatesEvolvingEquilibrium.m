% plot evolving equilibrium
clear
close all

nint = 10;  % number of interpolate states between two given original states
path = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/Evolving_equilibrium_35132els_10ts_P3/Lim2Div/';
load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/WEST_far_465.mat
load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/WEST_wall.mat
load([path,'Coordinates.mat'])
load([path,'Elements.mat'])
load([path,'fpsi.mat'])
load([path,'psir.mat'])
load([path,'psiz.mat'])
load([path,'psi.mat'])


nt = size(fpsi,2);
np = size(fpsi,1);
refEl = createReferenceElement(1,size(T,2));


% initialization
psi_int = zeros(np,(nt-1)*nint); 
fpsi_int = psi_int;
psir_int = psi_int;
psiz_int = psi_int;

xsi = linspace(0,1-1/nint,nint);


for i = 1:nt-1
    ind = (i-1)*nint + (1:nint);
    psi_int(:,ind)  = psi(:,i)*(1-xsi) + psi(:,i+1)*xsi;
    fpsi_int(:,ind) = fpsi(:,i)*(1-xsi) + fpsi(:,i+1)*xsi;
    psir_int(:,ind) = psir(:,i)*(1-xsi) + psir(:,i+1)*xsi;
    psiz_int(:,ind) = psiz(:,i)*(1-xsi) + psiz(:,i+1)*xsi;
end
psi = psi_int;
fpsi = fpsi_int;
psir = psir_int;
psiz = psiz_int;
save([path, 'fpsi_int_',num2str(nint),'.mat'],'fpsi')
save([path, 'psi_int_',num2str(nint),'.mat'],'psi')
save([path, 'psir_int_',num2str(nint),'.mat'],'psir')
save([path, 'psiz_int_',num2str(nint),'.mat'],'psiz')