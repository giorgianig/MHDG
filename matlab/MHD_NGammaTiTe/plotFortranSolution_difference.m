% plot Fortran solution
clear all
close all
global Mesh Mref
Mref = 12;

neq = 4;
nref = 5;
np = 10000;
lscale = 1.901*1e-3; 

% Circ
% solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
% meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
% solname = 'Sol_CircLimAlign_Quads_Nel208_P4_DPe0.100E+02_DPai0.100E+05_DPae0.100E+05_UNCP.h5';
% elemType = 0;

% West
meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';
solname1 = 'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.750E+01_DPai0.300E+06_DPae0.100E+07_NR0005.h5';
solname2 = 'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.750E+01_DPai0.300E+06_DPae0.100E+07_NR0006.h5';

elemType = 1;

% start
Mesh.lscale = lscale;
pos = strfind(solname1,'_P');
for i=1:10
    if strcmp(solname1(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname1(5:pos) '.h5'];

HDF5load([meshpath,meshname])

HDF5load([solpath,solname1])
u1 = transpose(reshape(u,[neq,numel(u)/neq]));
HDF5load([solpath,solname2])
u2 = transpose(reshape(u,[neq,numel(u)/neq]));

u = u1-u2;

if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end
    
    
refEl = createReferenceElement(elemType,size(T,2));


% up = cons2phys(u);  
% figure, plotSolution(X/lscale,T,up(:,1),refEl,nref);axis off,title('Density')
% figure, plotSolution(X/lscale,T,up(:,2),refEl,nref);axis off,title('Parallel velocity')
% figure, plotSolution(X/lscale,T,up(:,3),refEl,nref);axis off,title('Ions total energy')
% figure, plotSolution(X/lscale,T,up(:,4),refEl,nref);axis off,title('Electrons total energy')
% figure, plotSolution(X/lscale,T,up(:,5),refEl,nref);axis off,title('Ions pressure')
% figure, plotSolution(X/lscale,T,up(:,6),refEl,nref);axis off,title('Electrons pressure')
% figure, plotSolution(X/lscale,T,up(:,7),refEl,nref);axis off,title('Ions temperature')
% figure, plotSolution(X/lscale,T,up(:,8),refEl,nref);axis off,title('Electrons temperature')
% % figure, plotSolution(X/lscale,T,up(:,9),refEl,nref);axis off,title('Sound speed')
% figure, plotSolution(X/lscale,T,up(:,10),refEl,nref);axis off,title('Mach')

figure, plotSolution(X/lscale,T,u(:,1),refEl,nref);axis off,title('U1')
figure, plotSolution(X/lscale,T,u(:,2),refEl,nref);axis off,title('U2')
figure, plotSolution(X/lscale,T,u(:,3),refEl,nref);axis off,title('U3')
figure, plotSolution(X/lscale,T,u(:,4),refEl,nref);axis off,title('U4')

