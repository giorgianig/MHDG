% Postprocess routine

clear
close all

global testcase Mesh

%% Physical parameters
lscale_old = 1.901*1e-3;
lscale = 5.112369262649808e-04;

kT = 25;
testcase.n = 50;
axisym = 1;
fref = 7.28e6; % reference frequency
D = 0.02;
nref = 5;
npoints = 100;

%% Find the separatrix
load ../../Meshes/Meshes_2D/WEST/WEST_wall.mat
load ../../Meshes/Meshes_2D/WEST/WEST_far_465.mat
auxhugopar = load('/home/giorgio/Desktop/WEST_profs_s2d/parallel_profile.mat');
auxhugorad = load('/home/giorgio/Desktop/WEST_profs_s2d/radial_profile.mat');
auxMachp1 = load('/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/Mach_plus1.mat');
auxMachm1 = load('/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/Mach_moins1.mat');
auxMachHugo = ...
    [auxMachm1.c1;auxMachm1.c2;...
    auxMachm1.c3;auxMachm1.c4;...
    auxMachp1.c1;...
    auxMachp1.c2;auxMachp1.c3;...
    auxMachp1.c4;auxMachp1.c5;...
    auxMachp1.c6;auxMachp1.c7;...
    auxMachp1.c8;auxMachp1.c9];
%% Define separatrix
[Csep,h] = contour(r2D,z2D,flux2D,[-0.8662 -0.8662]);
close
aux = sqrt((Csep(1,:)-2.5).^2+Csep(2,:).^2)<sqrt((2.5-2.2722)^2+0.6082^2);
Csep = Csep(:,aux);
Centx = mean(Rwall);
Centy = mean(Zwall);
Csep_theta = pi/2+cart2pol(Csep(1,:)-Centx,Csep(2,:)-Centy);
Csep_theta(Csep_theta>pi)= Csep_theta(Csep_theta>pi)-2*pi;
[Csep_theta,ind] = sort(Csep_theta,'ascend');
Csep = Csep(:,ind);

%% Define another line
[Cline1,h] = contour(r2D,z2D,flux2D,[-0.86 -0.86]);
close
aux = sqrt((Cline1(1,:)-2.5).^2+(Cline1(2,:)+0.05).^2)< 0.75;
Cline1 = Cline1(:,aux);
df1 = diff(Cline1(1,:));
df2 = diff(Cline1(2,:));
dist = sqrt(df1.^2+df2.^2);
CurvCoo = [0,cumsum(dist)];
% Centx = mean(Rwall);
% Centy = mean(Zwall);
% Cline1_theta = pi/2+cart2pol(Cline1(1,:)-Centx,Cline1(2,:)-Centy);
% Cline1_theta(Cline1_theta>pi)= Cline1_theta(Cline1_theta>pi)-2*pi;
% [Cline1_theta,ind] = sort(Cline1_theta,'ascend');
% Cline1 = Cline1(:,ind);

% Reorder Cline1
% Cline_reo = zeros(size(Cline1));
% Cline_reo(:,1) = Cline1(:,1);
% for i=2:size(Cline1,2)-1
%     p = Cline1(:,i-1);
%     aux = Cline1(:,i:end);
%     
%     [~,ind]=min(sqrt((aux(1,:)-p(1)).^2+(aux(2,:)-p(2)).^2));
%     Cline_reo(:,i) = Cline1(:,ind+i-1);
% end
% Cline_reo(:,end) = Cline1(:,end);
% Cline1=Cline_reo;

%% Load mesh
path2Mesh = '../../Meshes/Meshes_2D/WEST/';
meshName = 'West_Hugo_h0.02_refCorn0.001_refSep0.01_P8.mat';
load([path2Mesh meshName])
refEl = createReferenceElement(1,size(T,2),[]);
X = X/lscale;
Mesh.X = X;
Mesh.T = T;
Mesh.lscale = lscale;

% Plot mesh with separatrix and line1
% plotMesh(X*lscale,T)
% plot(Rwall,Zwall,'k-'); axis equal
% hold on
% plot(Csep(1,:),Csep(2,:),'ko','Markersize',1)
% plot(Cline1(1,:),Cline1(2,:),'ko','Markersize',1)
% axis off
% readyforprintnew([8 6],24,[],[],1,[],[],'West_iso0.86')

%  stop

%% Load solution
path2Sol = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/';
solName = 'u_West_Hugo_h0.02_refCorn0.001_refSep0.01_P8_Diff.38000E-01.mat';
load([path2Sol solName])
% u = u0;
% % 
% % figure,plotSolution(X,T,u(1:2:end-1),refEl,5),axis off,title('Density')
% % readyforprintjpeg([8 6],24,[],[],[],[],[],'West_Dens_text')
% % 
% % figure,plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl,5),axis off,title('Mach')
% % readyforprintjpeg([8 6],24,[],[],[],[],[],'West_Mach_text')
% 
% figure,plotSolution(X,T,u(2:2:end)./u(1:2:end-1)/sqrt(kT),refEl,5,0,0,1),axis off
% hold on
% plot(auxMachHugo(:,1),auxMachHugo(:,2),'k*'),colorbar off
% axis([2.0581, 2.3969, -0.7824, -0.5153])
% p1 = plot(1,1,'ko');
% p2 = plot(1,1,'k*');
% legend([p1 p2],'HDG','SOLEDGE2D','location','south')
% readyforprintjpeg([8 6],24,[],[],[],[],[],'/home/giorgio/Dropbox/Conferences/PET_2017/figures/','MachZoneCompar_MHDGvsSOLEDGE')

% text(2.1,2.5,'HDG *')
% text(2.1,2.5,'SOLEDGE *')
% stop

%% Extract radial profile
% yhole = 10.^linspace(-12,0,100);
% xholel = 1162.27*lscale_old/lscale*ones(numel(yhole),1);
% xholer = 1518.19*lscale_old/lscale*ones(numel(yhole),1);
% line = [linspace(1.8,3,npoints)'/lscale,zeros(npoints,1)];
% uplot = evalDGapproximationAtPoints(line,u(1:2:end-1),X,T,refEl);
% uplot(uplot==0) = NaN;
% uplot(all([line(:,1)>xholel,line(:,1)<xholer],2)) =NaN;
% figure, semilogy(line(:,1),uplot,'k-','linewidth',2),grid on
% hold on
% auxhugorad.nmid(all([auxhugorad.rmid'*lscale_old/lscale>xholel(1),auxhugorad.rmid'*lscale_old/lscale<xholer(1)],2))=NaN;
% semilogy(auxhugorad.rmid*lscale_old/lscale,auxhugorad.nmid,'k--','linewidth',2)
% xlabel('$r$','interpreter','latex','fontname','times new roman','fontsize',20)
% ylabel('$n$','interpreter','latex','fontname','times new roman','fontsize',20)
% xlim([3500 5800])
% ylim([1e-12 1])
% hold on, semilogy(xholel,yhole,'k.')
% hold on, semilogy(xholer,yhole,'k.')
% leg = legend('HDG','SOLEDGE2D','Location','southeast');
% set(leg,'fontname','times new roman')
% readyforprintnew([8 6],24,[],[],1,[],[],'rho_West_Compar')
% readyforprintjpeg([8 6],24,[],[],[],[],[],'/home/giorgio/Dropbox/Conferences/PET_2017/figures/','RhoCurvesCompar')


%% Extract profile on separatrix
% line = Csep'/lscale;
% uplot = evalDGapproximationAtPoints(line,u(2:2:end),X,T,refEl);
% uplot(uplot==0) = NaN;
% figure, plot3(line(:,1)*lscale,line(:,2)*lscale,uplot,'r-'),grid on

%% Extract profile on separatrix
line = Cline1'/lscale;
uplot = evalDGapproximationAtPoints(line,u(2:2:end)./u(1:2:end-1)/sqrt(kT),X,T,refEl);
uplot(uplot==0) = NaN;
% figure, plot3(line(:,1)*lscale,line(:,2)*lscale,uplot,'ro'),grid on
CurvCoo = CurvCoo(~isnan(uplot));
uplot = uplot(~isnan(uplot));
CurvCoo = CurvCoo-min(CurvCoo);
figure, plot(CurvCoo/lscale,uplot,'k-','linewidth',2),grid on
hold on
plot(auxhugopar.s*lscale_old/lscale,auxhugopar.M,'k--','linewidth',2)
legend('HDG','SOLEDGE2D','Location','southeast')

xlabel('$s$','interpreter','latex','fontname','times new roman','fontsize',20)
ylabel('$M_{\|}$','interpreter','latex','fontname','times new roman','fontsize',20)
xlim([0 6600])
% readyforprintnew([8 6],24,[],[],1,[],[],'upar_West_Compar')
% readyforprintjpeg([8 6],24,[],[],[],[],[],'/home/giorgio/Dropbox/Conferences/PET_2017/Presentation/','MachCurvesCompar')


