% plot Fortran solution
clear
close all
global Mesh Mref neq theta ntor
Mref = 12;
%  clc
neq = 4;
nref = 10;
np = 10000;
lscale = 1.901*1e-3;

n0 = 1e19;
T0 = 50;


ntpos =2;
theta = 2;        % Toroidal angle considered

dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);


ntpos=1;
tpos = tpos(1);
tpos = 1;



% Circ -quads
% solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/Circ/';
% meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/Circ/';
% solname = 'Sol_CircLimAlign_Quads_Nel2916_P4_DPe0.100E+01_DPai0.300E+06_DPae0.100E+08.h5';
% elemType = 0;

% Circ -tria
% solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
% meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
% solname = 'Sol_Circle_LIM_InfThin_h0.15_RefLim0.02_P8_DPe0.380E-01_DPai0.300E+06_DPae0.100E+08.h5';
% elemType = 1;

% West
% meshpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
% solpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
% solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P4_DPe0.240E+00_DPai0.300E+06_DPae0.100E+08.h5';
% solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_NoHole_P3_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08_NR0005.h5';
% solname = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P7_DPe0.398E-01_DPai0.300E+06_DPae0.100E+08_NR0050.h5';
% solname = 'Sol2D_CircLimAlign_Quads_Nel208_P4_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08.h5';
% solname = 'Sol2D_Circle_LIM_InfThin_h0.1_RefCorn0.02_P4_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08_NR0039.h5';
% solname = 'Sol2D_Circle_lim_h02_P4_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08_NR0025.h5';
% elemType = 1;

% Square
meshpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
% solname = 'Sol2D_Square_3_P3_DPe0.100E+01';
% solname = 'Sol3D_Square_3_P3_Ntor1Ptor1_DPe0.100E+01.h5';
solname = 'Sol2D_CircLimAlign_Quads_Nel208_P4_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08_NR0001';
% solname = 'Sol3D_Square_3_P3_Ntor1Ptor1_DPe0.100E+01_DPai0.200E+02_DPae0.300E+02';
% elemType = 1;
% lscale = 1;


% solname = 'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.380E-01_DPai0.300E+06_DPae0.100E+08_0040.h5';
% solname = 'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.100E+00_DPai0.300E+06_DPae0.100E+07_NR0011.h5';
% solname = 'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.250E+00_DPai0.300E+06_DPae0.100E+07_NR0038.h5';
% solname =  'Sol_West_Hugo_h0.04_refCorn_P3_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08.h5';

% start
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
u = transpose(reshape(u,[neq,numel(u)/neq]));
if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end

refEl = createReferenceElement(elemType,size(T,2));
% plotMesh(X,T,elemType),axis off


%% load magnetic field
% magneticName = [meshpath 'MagneticField'  meshname(1:end-3) '.mat'];
% load(magneticName)
for ieq = 1:neq
    
    if strcmpi(solname(4:5),'2D')
        figure()
        if ieq~=2
            if ieq==1
                plotSolution(X/lscale,T,u(:,ieq),refEl,5)
            else
                plotSolution(X/lscale,T,u(:,ieq)./u(:,1),refEl,5)
            end
        else
            plotSolution(X/lscale,T,u(:,ieq)./u(:,1),refEl,5)
        end
        hold on, plotMesh(X/lscale,T,elemType,[])
        title(['U',num2str(ieq)])
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
            figure(),clf
            upol = extractSolutionInAtGivenTheta(u(:,1),T,refEl,refElTor,tpos(itor));
            plotSolution(X/lscale,T,upol,refEl)
            %     caxissave(itor,:) = caxis;
            %         title(['Plane ' num2str(itor)])
            drawnow
        end
        
        
    end
end





% hold on, plot(X(Tb(boundaryFlag==7,:),1),X(Tb(boundaryFlag==7,:),2),'r-','linewidth',1)
% print('-djpeg','-r600','/home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MeshTria.jpg')
% system(['convert /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MeshTria.jpg -trim /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MeshTria.jpg'])
% stop

% figure, plotSolution(X/lscale,T,col(transpose(scdiff_nodes)),refEl,nref);axis off,title('Shock-capturing diffusion')
% stop


% u = filterSolution(T,u,refEl);








% up = cons2phys(u);
% figure, plotSolution(X/lscale,T,up(:,1)*n0,refEl,nref);axis off,title('Density')
% % hLegend = findobj(gcf, 'Type', 'colorbar');
% % set(hLegend, 'fontsize',16,'fontName','times new roman')
% % print('-djpeg','-r600','/home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/DensWest.jpg')
% % system(['convert /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/DensWest.jpg -trim /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/DensWestRef.jpg'])
% % stop
% % % plotMesh(X/lscale,T,1),axis(1.0e+03 *[1.2758    1.2878    0.4071    0.4166])
% % % caxis([0.0108    1.0001])
% figure, plotSolution(X/lscale,T,up(:,2),refEl,nref);axis off,title('Parallel velocity')
% figure, plotSolution(X/lscale,T,up(:,3),refEl,nref);axis off,title('Ions total energy')
% figure, plotSolution(X/lscale,T,up(:,4),refEl,nref);axis off,title('Electrons total energy')
% figure, plotSolution(X/lscale,T,up(:,5),refEl,nref);axis off,title('Ions pressure')
% figure, plotSolution(X/lscale,T,up(:,6),refEl,nref);axis off,title('Electrons pressure')
% figure, plotSolution(X/lscale,T,up(:,7)*T0,refEl,nref);axis off,title('Ions temperature'),
% % hLegend = findobj(gcf, 'Type', 'colorbar');
% % set(hLegend, 'fontsize',16,'fontName','times new roman')
% % print('-djpeg','-r600','/home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/IonTempWest.jpg')
% % system(['convert /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/IonTempWest.jpg -trim /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/IonTempWestRef.jpg'])
% % stop
%
% % caxis([0.1658    1.0001])
% % hold on, plotMesh(X/lscale,T,1),axis(1.0e+03 *[1.2758    1.2878    0.4071    0.4166])
% figure, plotSolution(X/lscale,T,up(:,8)*T0,refEl,nref);axis off,title('Electrons temperature')
% % hLegend = findobj(gcf, 'Type', 'colorbar');
% % set(hLegend, 'fontsize',16,'fontName','times new roman')
% % print('-djpeg','-r600','/home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/EleTempWest.jpg')
% % system(['convert /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/EleTempWest.jpg -trim /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/EleTempWestRef.jpg'])
% % stop
% figure, plotSolution(X/lscale,T,up(:,9),refEl,nref);axis off,title('Sound speed')
% figure, plotSolution(X/lscale,T,up(:,10),refEl,nref);axis off,%title('Mach')
% hLegend = findobj(gcf, 'Type', 'colorbar');
% set(hLegend, 'fontsize',16,'fontName','times new roman')
% print('-djpeg','-r600','/home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MachWest.jpg')
% system(['convert /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MachWest.jpg -trim /home/giorgio/Dropbox/Papers/Paper_Miei/Paper_PSI/figures/newFig/MachWestRef.jpg'])
% stop

% figure, plotSolution(X/lscale,T,u(:,1),refEl,nref);axis off,title('U1')
% figure, plotSolution(X/lscale,T,u(:,2),refEl,nref);axis off,title('U2')
% figure, plotSolution(X/lscale,T,u(:,3),refEl,nref);axis off,title('U3'),
% caxis([0.0324   18.0014])
% hold on, plotMesh(X/lscale,T,1),axis(1.0e+03 *[1.2758    1.2878    0.4071    0.4166])
% hold on, plot(X(T(1902,:),1)/lscale,X(T(1902,:),2)/lscale,'r*')
% hold on, text(X(T(1902,:),1)/lscale,X(T(1902,:),2)/lscale,num2str(transpose(1:size(T,2))),'fontsize',16)
% figure, plotSolution(X/lscale,T,u(:,4),refEl,nref);axis off,title('U4')

% figure, plotMesh(X,T,1),hold on
% velx = u(:,2).*col(Magnetic.bxnodes);
% vely = u(:,2).*col(Magnetic.bynodes);
% x = X(:,1);
% y = X(:,2);
% x = col(x(T'));
% y = col(y(T'));
% quiver(x,y,velx,vely,1)

% plotSolutionPhys(X/lscale,T,u,refEl,nref)

% rho = 0.9;
% theta = linspace(0,2*pi,np)+pi/2;
% x = rho*cos(theta);
% y = rho*sin(theta);
% line = [x',y'];
% uplot = evalDGapproximationAtPoints(line,u,X,T,refEl);
% up_plot = cons2phys(uplot);
% figure,plot(theta-3/2*pi,up_plot(:,7)),title('Ions temperature')
% figure,plot(theta-3/2*pi,up_plot(:,8)),title('Electrons temperature')


% iel = 2196;
% nv = size(refEl.NodesCoord,1);
% ind = (iel-1)*nv+(1:nv);
% disp('pressure')
% up(ind,5)
% disp('density')
% up(ind,1)
% disp('temperature')
% up(ind,5)./up(ind,1)
%
% iel = 2210;
% nv = size(refEl.NodesCoord,1);
% ind = (iel-1)*nv+(1:nv);
% disp('pressure')
% up(ind,5)
% disp('density')
% up(ind,1)
% disp('temperature')
% up(ind,5)./up(ind,1)


% iel = 31685;
% nv = size(refEl.NodesCoord,1);
% ind = (iel-1)*nv+(1:nv);
% disp('density')
% up(ind,1)
% disp('pressure')
% up(ind,5)
% disp('temperature')
% up(ind,7)
% disp('U3')
% u(ind,3)
%
% disp('U3')
% u(ind,3)
% disp('pressure')
% up(ind,5)
% disp('density')
% up(ind,1)
% disp('temperature')
% up(ind,5)./up(ind,1)
%
%
%
% intsol = computeIntSol(X,T,u,refEl);
% disp(intsol(1902,3))