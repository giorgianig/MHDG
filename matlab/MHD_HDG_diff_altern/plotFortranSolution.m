% plot Fortran solution
clear 
% close all
global Mesh

global   theta 
global    ntor ptor refElPol refElTor 


kT = 25;
lscale = 1.901*1e-3; 
% solpath = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/NoHole/';
% solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/';
meshpath = solpath;% '/home/giorgio/Desktop/MHDG_Thomas/test/Evolv_Giorgio/';
% solpath = '/home/giorgio/Saves_MHDG_Marconi/Goldston/WithDrift/';
% solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P8_Diff.38000E-01_NoDrift_200eV.h5';
% solname = 'Sol_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.003_refSep0.008_P10_Diff.31623E-04_0006.h5';
% solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P8_Diff.38000E-01.h5';
solname = 'Sol3D_CircLimAlign_Quads_Nel208_P4_Ntor2Ptor1_DPe0.380E+00.h5';
% solname = 'Sol_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P5_Diff.23714E-03_0034.h5';
% start


solname = 'Sol3D_Square_3_P4_Ntor2Ptor1_DPe0.300E+01.h5';

ntpos =2;
theta = 2*pi;        % Toroidal angle considered
dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);



Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(7:pos) '.h5'];
HDF5load([meshpath,meshname])
HDF5load([solpath,solname])
u = u';

if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end

refEl = createReferenceElement(elemType,size(T,2));



if strcmpi(solname(4:5),'2D')
    figure, plotSolution(X/lscale,T,u(1:2:end-1),refEl,5);axis off
    figure, plotSolution(X/lscale,T,u(2:2:end)./u(1:2:end-1)./sqrt(kT),refEl,5);axis off
    hold on, plotMesh(X,T,elemType,[])
elseif strcmpi(solname(4:5),'3D')
    k = strfind(solname,'Ntor');
    for ipos = 1:3
       if strcmpi(solname(k+4+ipos),'P'),break,end  
    end
    ntor = str2double(solname(k+4:k+4+ipos-1));
    k = strfind(solname,'Ptor');
    ptor = str2double(solname(k+4));
    refElTor = createReferenceElement(0,(ptor+1)^2,[]);
    for itor = 1:ntpos
        upol1 = extractSolutionInAtGivenTheta(u(1:2:end-1),T,refEl,refElTor,tpos(itor));
        figure(),clf
        plotSolution(X/lscale,T,upol1,refEl);axis off
        title(['U1 - Plane ' num2str(itor)])
        drawnow
        upol2 = extractSolutionInAtGivenTheta(u(2:2:end),T,refEl,refElTor,tpos(itor));
        figure(),clf
%         plotSolution(X/lscale,T,upol2./upol1./sqrt(kT),refEl);axis off
        plotSolution(X/lscale,T,upol2,refEl);axis off
        title(['U2 - Plane ' num2str(itor)])
        drawnow        
    end
end



