% plot Fortran solution
clear all
close all
global Mesh Mref
Mref = 12;
T0 =50;                 % Temperature (eV)
% mi = 3.35e-27;       % Ionic mass (kg)
% me = 9.109e-31;    % electronic mass (kg)
% B0 = 1;                   % Magnetic field at magnetic axis (Tesla)
% D = 300;                 % Diffusion in m^2/s
% k0 = 2000;             % Stangeby
n0 = 1e19;             % Reference density

neq = 4;
nref = 5;

lscale = 1.901*1e-3;
solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';

np = 1000;
%% Solution 1
meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
elemType = 1; % 1-Triangles/0-Quadrangles
solname = 'Sol_Circle_LIM_InfThin_h0.15_RefLim0.02_P8_DPe0.380E-01_DPai0.300E+06_DPae0.100E+08.h5';


% start
Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];

HDF5load([meshpath,meshname])
HDF5load([solpath,solname])
u1 = transpose(reshape(u,[neq,numel(u)/neq]));
refEl = createReferenceElement(1,size(T,2));
up1 = cons2phys(u1);
figure(1), plotMesh(X,T,1), axis off
print('-djpeg','-r600','/home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshTria.jpg')
system(['convert /home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshTria.jpg -trim /home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshTria.jpg'])


linex = linspace(0.5,1,np)';
linex = linex/lscale;
liney = zeros(size(linex));
line = [linex,liney];
uplot = evalDGapproximationAtPoints(line,n0*up1(:,1),X/lscale,T,refEl); % density
figure(3), semilogy(linex,uplot,'b'),grid on
uplot = evalDGapproximationAtPoints(line,T0*up1(:,7),X/lscale,T,refEl); % ions temperature
figure(4), semilogy(linex,uplot,'b'),grid on
uplot = evalDGapproximationAtPoints(line,T0*up1(:,8),X/lscale,T,refEl); % electrons temperature
figure(5), semilogy(linex,uplot,'b'),grid on



%% Solution 2
meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/';
elemType = 0; % 1-Triangles/0-Quadrangles
solname = 'Sol_CircLimAlign_Quads_Nel208_P4_DPe0.380E-01_DPai0.300E+06_DPae0.100E+08.h5';


% start
Mesh.lscale = lscale;
pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];

HDF5load([meshpath,meshname])
HDF5load([solpath,solname])
u2 = transpose(reshape(u,[neq,numel(u)/neq]));
refEl = createReferenceElement(0,size(T,2));
figure(2), plotMesh(X,T,0), axis off
print('-djpeg','-r600','/home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshQua.jpg')
% system(['convert ' path name '.jpg -trim ' path  name '.jpg'])
system(['convert /home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshQua.jpg -trim /home/giorgio/Documents/EoCoE/WP1_applied_research/tech_report/HOFem_Tokam3X/figures/MeshQua.jpg'])
up2 = cons2phys(u2);

linex = linspace(0.5,1,np)';
linex = linex/lscale;
liney = zeros(size(linex));
line = [linex,liney];

uplot = evalDGapproximationAtPoints(line,n0*up2(:,1),X/lscale,T,refEl); % density
figure(3), hold on, semilogy(linex,uplot,'r'),grid on
legend('Tria','Quads','location','southeast','fontname','times new roman')
xlabel('$r$','fontname','times new roman','interpreter','latex')
ylabel('$n [m^{-3}]$ ','fontname','times new roman','interpreter','latex')
set(gca,'fontname','times new roman')
readyforprintnew([8 6],24,[],[],1,[],[],'density')

uplot = evalDGapproximationAtPoints(line,T0*up2(:,7),X/lscale,T,refEl); % ions temperature
figure(4), hold on, semilogy(linex,uplot,'r'),grid on
legend('Tria','Quads','location','southeast','fontname','times new roman')
xlabel('$r$','fontname','times new roman','interpreter','latex')
ylabel('Ions temperature [eV]','fontname','times new roman','interpreter','latex')
set(gca,'fontname','times new roman')
readyforprintnew([8 6],24,[],[],1,[],[],'ionsTemp')

uplot = evalDGapproximationAtPoints(line,T0*up2(:,8),X/lscale,T,refEl); % electrons temperature
figure(5), hold on, semilogy(linex,uplot,'r'),grid on
legend('Tria','Quads','location','southeast','fontname','times new roman')
xlabel('$r$','fontname','times new roman','interpreter','latex')
ylabel('Electrons temperature [eV]','fontname','times new roman','interpreter','latex')
set(gca,'fontname','times new roman')
 readyforprintnew([8 6],24,[],[],1,[],[],'electTemp')
