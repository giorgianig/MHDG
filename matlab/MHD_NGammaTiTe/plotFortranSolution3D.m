% plot Fortran solution
clear
close all
global Mesh Mref neq theta ntor
Mref = 12.5;
%  clc
neq = 2; % number of equations solved
nref = 10;
np = 10000;
lscale = 1.901*1e-3;

n0 = 1e19;
T0 = 50;

% Simulations parameters
tmax = 2*pi/9; %2*pi/9.;
R0 = 3.4; %2.5; %3.4;

% 3D only
ntpos =17; % number of toroidal position to plot
theta = tmax; % Toroidal angle considered
dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);
row_subplot = 2;
col_subplot = 2;
consSol = false; % set to true if plot for conserved fields
mesh2plot = false; %set to true to plot mesh on the solution
fields2plot = [1,2]; % if consSol 1: N, 2: Gamma , 3: NEi, 4: NEe ELSE 1: N, 2: u para (aniso)/ Mach(iso), 3: Ei, 4: Ee, 5: Pi, 6: Pe, 7: Ti, 8: Te, 9: Cs, 10: Mach
B3D = true; % set to true if 3D magnetic field
diff = true; % true to compute difference with ref
average = false; % true to compute difference with toroidal average
str_diff = 'diff_ref'; % diff_ref or diff_av or nodiff for saving
saveplot = false;

fig = 1;

% Circ -quads
% elemType = 0;

% Circ -tria
% elemType = 1;

%home = '/home/bluce/MHDG_sim/iso/3D/';
home = '/home/laurence/Documents/Benjamin/MHDG/iso/3D/';
%home = '/home/bluce/Documents/phd-dev/MHDG/test/test_dir/';
%savedir = '/home/bluce/Images/Octave-Matlab/MHDG/';
savedir = '/home/laurence/Documents/Benjamin/Images/MHDG/';

%home ='/Home/BL254371/Simulations/MHDG/';
%savedir = '/Home/BL254371/Images/MHDG/';

% ref iso square
% ref = 'InfLim909_P5/ref/rmp/';
ref = 'InfLim909_P5/ref/ripple/';
% ref = 'OriLim843_P5/ref/rmp/';
% ref = 'OriLim843_P5/ref/ripple/';
% ref = 'VerLim907_P5/ref/rmp/';
% ref = 'VerLim907_P5/ref/ripple/';
% ref = 'West3043_P5/ref/rmp/';
% ref = 'West3043_P5/ref/ripple/';
meshpath = [home,ref];
solpath = [home,ref];
coilspath = [home,ref];
% solname = 'Sol3D_Circ_InfLIM_Quads_YesHole_Nel909_P5_Ntor3Ptor4_DPe0.100E+01.h5';
solname = 'Sol3D_Circ_InfLIM_Quads_YesHole_Nel909_P5_Ntor2Ptor4_DPe0.100E+01.h5';
% solname = 'Sol3D_Circ_OriLIM_Quads_YesHole_Nel843_P5_Ntor3Ptor4_DPe0.380E-01.h5';
% solname = 'Sol3D_Circ_OriLIM_Quads_YesHole_Nel843_P5_Ntor2Ptor4_DPe0.380E-01.h5';
% solname = 'Sol3D_Circ_VerLIM_Quads_YesHole_Nel907_P5_Ntor3Ptor4_DPe0.380E-01.h5';
% solname = 'Sol3D_Circ_VerLIM_Quads_YesHole_Nel907_P5_Ntor2Ptor4_DPe0.380E-01.h5';
% solname = '';
% solname = 'Sol3D_West_YesHole_Nel3043_P5_Ntor2Ptor4_DPe0.100E+00.h5';


coilsname = '';
if diff
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
    % Load mesh and solution
    HDF5load([meshpath,meshname])
    HDF5load([solpath,solname])
    % u contains N, Gamma, Ti and Te (from u(:,1) to u(:,4))
    u0 = transpose(reshape(u,[neq,numel(u)/neq]));
    up0 = cons2phys(u0);
end

% % rmp iso square pi/2
% nbCoils = 8;
% nbRow = 2;
% rmp = 'InfLim909_P5/rmp/';
% % rmp = 'OriLim843_P5/rmp/';
% % rmp = 'VerLim907_P5/rmp/';
% % rmp = 'West3043_P5/rmp/';
% meshpath = [home,rmp];
% solpath = [home,rmp];
% coilspath = [home,rmp];
% solname = 'Sol3D_Circ_InfLIM_Quads_YesHole_Nel909_P5_Ntor3Ptor4_DPe0.100E+01.h5';
% % solname = 'Sol3D_Circ_OriLIM_Quads_YesHole_Nel843_P5_Ntor3Ptor4_DPe0.380E-01.h5';
% % solname = 'Sol3D_Circ_VerLIM_Quads_YesHole_Nel907_P5_Ntor3Ptor4_DPe0.380E-01.h5';
% % solname = '';
% coilsname = {};
% for i=1:nbRow
%     str = ['RMP_',pad(int2str(nbCoils),2,'left','0'),'_coils_',pad(int2str(i),2,'left','0'),'_',pad(int2str(nbRow),2,'left','0'),'.h5'];
%     coilsname = [coilsname,str];
% end
% savedir = [savedir,rmp];
% %mkdir(savedir);

% ripple iso square pi/8 (16 coils) or 2*pi/9 (18 coils)
nbCoils = 18;
nbRow = 1;
ripple = 'InfLim909_P5/ripple/';
% ripple = 'OriLim843_P5/ripple/';
% ripple = 'VerLim907_P5/ripple/';
% ripple = 'West3043_P5/ripple/';
meshpath = [home,ripple];
solpath = [home,ripple];
coilspath = [home,ripple];
solname = 'Sol3D_Circ_InfLIM_Quads_YesHole_Nel909_P5_Ntor2Ptor4_DPe0.100E+01.h5';
% solname = 'Sol3D_Circ_OriLIM_Quads_YesHole_Nel843_P5_Ntor2Ptor4_DPe0.380E-01.h5';
% solname = 'Sol3D_Circ_VerLIM_Quads_YesHole_Nel907_P5_Ntor2Ptor4_DPe0.380E-01.h5';
% solname = 'Sol3D_West_YesHole_Nel3043_P5_Ntor2Ptor4_DPe0.100E+00.h5';
coilsname = {};
for i=1:nbRow
    str = ['Ripple_',pad(int2str(nbCoils),2,'left','0'),'_coils.h5'];
    coilsname = [coilsname,str];
end
savedir = [savedir,ripple];
%mkdir(savedir);

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

% Load mesh and solution
HDF5load([meshpath,meshname])
HDF5load([solpath,solname])
% u contains N, Gamma, Ti and Te (from u(:,1) to u(:,4))
u = transpose(reshape(u,[neq,numel(u)/neq]));
up = cons2phys(u);

% Some calc for size of the simulation box
s_x = size(X,1);
xlim1 = min(X(:,1)/lscale) - 0.05*(max(X(:,1)/lscale) - min(X(:,1)/lscale));
xlim2 = max(X(:,1)/lscale) + 0.05*(max(X(:,1)/lscale) - min(X(:,1)/lscale));
ylim1 = min(X(:,2)/lscale) - 0.05*(max(X(:,2)/lscale) - min(X(:,2)/lscale));
ylim2 = max(X(:,2)/lscale) + 0.05*(max(X(:,2)/lscale) - min(X(:,2)/lscale));

if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end

refEl = createReferenceElement(elemType,size(T,2));

if average
    k = strfind(solname,'Ntor');
    for ipos = 1:3
        if strcmpi(solname(k+4+ipos),'P'),break,end
    end
    ntor = str2double(solname(k+4:k+4+ipos-1));
    k = strfind(solname,'Ptor');
    ptor = str2double(solname(k+4));
    refElTor = createReferenceElement(0,(ptor+1)^2,[]);
    sz = size(extractSolutionInAtGivenTheta(u(:,1),T,refEl,refElTor,tpos(1)),1);
    um = zeros(ntpos,sz,size(u,2));
    upm = zeros(ntpos,sz,size(up,2));
    for ieq = fields2plot
        for itor = 1:ntpos
            um(itor,:,ieq) = extractSolutionInAtGivenTheta(u(:,ieq),T,refEl,refElTor,tpos(itor));
            upm(itor,:,ieq) = extractSolutionInAtGivenTheta(up(:,ieq),T,refEl,refElTor,tpos(itor));
        end
    end
    %stop
    um = squeeze(mean(um,1));
    upm = squeeze(mean(upm,1));
    %stop
end
% if diff
%     %Comparison
%     % u = (u-u0);
%     % up = (up-up0);
%     u(:,1) = (u(:,1)-u0(:,1))*100./u0(:,1);
%     up(:,1) = (up(:,1)-up0(:,1))*100./up0(:,1);
%     u(:,2) = (u(:,2)-u0(:,2));
%     up(:,2) = (up(:,2)-up0(:,2));
%     test = abs(u) > 500;
%     u(test) = nan;
%     test = abs(up) > 500;
%     up(test) = nan;
% end


%% Mesh
if mesh2plot
    figure(fig), fig = fig+1;
    plotMesh(X,T,elemType)%,'plotNodes')
    axis on
    box on
    if saveplot
        saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_Mesh.png'])
    end
end

% Mesh 3D
% Load magnetic coils coordinates and plot
if B3D
    % Load magnetic field
    HDF5load([coilspath,'Bfield.h5'])
    BR = mat(:,1);
    BZ = mat(:,2);
    BPhi = mat(:,3);
    Bmod = sqrt(mat(:,1).^2 + mat(:,2).^2 + mat(:,3).^2);
    HDF5load([coilspath,'Bperturbated.h5'])
    BRpert = mat(:,1);
    BZpert = mat(:,2);
    BPhipert = mat(:,3);
    Bmodpert = sqrt(mat(:,1).^2 + mat(:,2).^2 + mat(:,3).^2);
    
    if elemType == 0
        faceNodes = faceNodes_aux_quads(size(T,2));
    elseif elemType == 1
        faceNodes = faceNodes_aux(size(T,2));
    end
    optn = 1;
    %Ordering the face nodes in a row vector without connectivity between them
    [nOfFaces,nOfNodesPerFace] = size(faceNodes);
    oFaceNodes = zeros(1,nOfFaces*(nOfNodesPerFace-1));
    np = nOfNodesPerFace - 1;
    aux = 1 - np;
    aux2 = 0;
    for iface = 1:nOfFaces
        aux = aux + np;
        aux2 = aux2 + np;
        oFaceNodes(aux:aux2) = faceNodes(iface,1:np);
    end
    nt = 32; % toroidal number of mesh to plot
    smat = size(X,1)*nt;
    X3D = zeros(smat,3);
    T3D = zeros(size(T,1)*nt,size(T,2));
    
    r = sqrt((X(:,1) - R0).^2 + X(:,2).^2);
    Z = X(:,2);
    R = r.*cos(asin(Z./r));
    
    for i=1:nt
        tor_angle = tmax*i/nt;
        X3D((i-1)*size(X,1)+1:i*size(X,1),1) = R.*cos(tor_angle);
        X3D((i-1)*size(X,1)+1:i*size(X,1),2) = R.*sin(tor_angle);
        X3D((i-1)*size(X,1)+1:i*size(X,1),3) = Z;
        T3D((i-1)*size(T,1)+1:i*size(T,1),:) = T + (i-1)*size(X,1);
    end
    %Connectivity for the faces
    patchFaces = T3D(:,oFaceNodes);
    %Plot 3D mesh
    figure(fig), fig = fig+1;
    patchHandle = patch('Faces',patchFaces,'Vertices',X3D,'FaceColor','none','EdgeAlpha',1,'linewidth',0.4, 'linestyle',':');
    
    % Plot of the coils
    % Load coils
    for l=1:nbRow
        HDF5load([coilspath,cell2mat(coilsname(l))]);
        n = numel(mat(:,1)); % 2 positions per line (3 coordinates)
        mat = mat.*lscale;
        for i=1:n % segments
            line([mat(i,1) mat(i,2)],[mat(i,3) mat(i,4)],[mat(i,5) mat(i,6)], 'color', 'red', 'linewidth',2);
        end
    end 
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    box on;
%     if saveplot
%         saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_Coils.png'])
%     end
    %stop
    % plot pertubated magnetic field
    % Module
    %figure(fig), fig = fig+1;
    clim1 = min(Bmodpert);
    clim2 = max(Bmodpert);
    j = 1;
    step = 1;
    tposB = linspace(0+dth,theta-dth,size(Bmodpert,1)/s_x/step);
    for i = 1:step:size(Bmodpert,1)/s_x
        figure(fig), fig = fig+1;
        %subplot(3,5,j)
        %tiledlayout(3,5)
        %nexttile
        title(['|Bpert|, plane :',num2str(tposB(j)),'/',num2str(tmax)])
        plotSolution(X/lscale,T,Bmodpert((i-1)*s_x+1:i*s_x),refEl,5);
        xlim([xlim1,xlim2]);
        ylim([ylim1,ylim2]);
        caxis([clim1, clim2]);
        box on;
        if saveplot
            saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_Bert_',pad(int2str(j),3,'left','0'),'.png'])
        end
        j = j + 1;
    end
    stop
    % BR
    %figure(fig), fig = fig+1;
    clim1 = min(BRpert);
    clim2 = max(BRpert);
    j = 1;
    step = 1;
    tposB = linspace(0+dth,theta-dth,size(BRpert,1)/s_x/step);
    for i = 1:step:size(BRpert,1)/s_x
        %subplot(3,5,j)
        figure(fig), fig = fig+1;
        %tiledlayout(3,5)
        %nexttile
        title(['|BRpert|, plane :',num2str(tposB(j)),'/',num2str(tmax)])
        plotSolution(X/lscale,T,BRpert((i-1)*s_x+1:i*s_x),refEl,5);
        xlim([xlim1,xlim2]);
        ylim([ylim1,ylim2]);
        caxis([clim1, clim2]);
        box on;
        if saveplot
            saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_BRert_',pad(int2str(j),3,'left','0'),'.png'])
        end
        j = j + 1;
    end
    % BZ
    %figure(fig), fig = fig+1;
    clim1 = min(BZpert);
    clim2 = max(BZpert);
    j = 1;
    step = 1;
    tposB = linspace(0+dth,theta-dth,size(BZpert,1)/s_x/step);
    for i = 1:step:size(BZpert,1)/s_x
        figure(fig), fig = fig+1;
        %subplot(3,5,j)
        %tiledlayout(3,5)
        %nexttile
        title(['|BZpert|, plane :',num2str(tposB(j)),'/',num2str(tmax)])
        plotSolution(X/lscale,T,BZpert((i-1)*s_x+1:i*s_x),refEl,5);
        xlim([xlim1,xlim2]);
        ylim([ylim1,ylim2]);
        caxis([clim1, clim2]);
        box on;
        if saveplot
            saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_BZert_',pad(int2str(j),3,'left','0'),'.png'])
        end
        j = j + 1;
    end
    % BPhi
    %figure(fig), fig = fig+1;
    clim1 = min(BPhipert);
    clim2 = max(BPhipert);
    j = 1;
    step = 1;
    tposB = linspace(0+dth,theta-dth,size(BPhipert,1)/s_x/step);
    for i = 1:step:size(BPhipert,1)/s_x
        figure(fig), fig = fig+1;
        %subplot(3,5,j)
        %tiledlayout(3,5)
        %nexttile
        title(['|BPhipert|, plane :',num2str(tposB(j)),'/',num2str(tmax)])
        plotSolution(X/lscale,T,BPhipert((i-1)*s_x+1:i*s_x),refEl,5);
        xlim([xlim1,xlim2]);
        ylim([ylim1,ylim2]);
        caxis([clim1, clim2]);
        box on;
        if saveplot
            saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_BPhiert_',pad(int2str(j),3,'left','0'),'.png'])
        end
        j = j + 1;
    end
end
%stop
%Plot solution
%if conSol 1: N, 2: Gamma, 3: NEi, 4: NEe
%else 1: N, 2: u para, 3: Ei, 4: Ee, 5: Pi, 6: Pe, 7: Ti, 8: Te
if strcmpi(solname(4:5),'3D') %3D
    close all
    if consSol
        k = strfind(solname,'Ntor');
        for ipos = 1:3
            if strcmpi(solname(k+4+ipos),'P'),break,end
        end
        ntor = str2double(solname(k+4:k+4+ipos-1));
        k = strfind(solname,'Ptor');
        ptor = str2double(solname(k+4));
        refElTor = createReferenceElement(0,(ptor+1)^2,[]);
        for ieq = fields2plot
            %figure(fig), fig = fig+1;
            %tiledlayout(row_subplot,col_subplot)
            clim1 = min(u(:,ieq));
            clim2 = max(u(:,ieq));
            for itor = 1:ntpos
                figure(fig), fig = fig+1;
                %subplot(row_subplot,col_subplot,itor)
                %nexttile
                upol = extractSolutionInAtGivenTheta(u(:,ieq),T,refEl,refElTor,tpos(itor));
                if diff
                    upol0 = extractSolutionInAtGivenTheta(u0(:,ieq),T,refEl,refElTor,tpos(itor));
                    upol0bis = extractSolutionInAtGivenTheta(u0(:,ieq),T,refEl,refElTor,tpos(itor));
                elseif average
                    upol0 = um(:,ieq);
                    upol0bis = um(:,ieq);
                else
                    upol0 = 0.;
                    upol0bis = 1.;
                end
                if ieq == 1
                    upol = (upol - upol0)*100./upol0bis;
                    title(['N, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'N';
                    %clim1 = -8;
                    %clim2 = 4;
                elseif ieq == 2
                    upol = (upol - upol0);
                    title(['Gamma, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Gamma';
                elseif ieq == 3
                    upol = (upol - upol0)*100./upol0bis;
                    title(['NEi, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'NEi';
                elseif ieq == 4
                    upol = (upol - upol0)*100./upol0bis;
                    title(['NEe, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'NEe';
                else
                    upol = (upol - upol0);
                    title(['???, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'bob';
                end               
                plotSolution(X/lscale,T,upol,refEl);
                if mesh2plot
                    hold on, plotMesh(X/lscale,T,elemType,[]);
                end
                xlim([xlim1,xlim2]);
                ylim([ylim1,ylim2]);
                caxis([clim1,clim2]);
                box on;
                %caxissave(itor,:) = caxis;
                drawnow
                if saveplot
                    saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_',str_diff,'_',f,'_',pad(int2str(itor),3,'left','0'),'.png'])
                end
            end
        end
    else
        close all
        k = strfind(solname,'Ntor');
        for ipos = 1:3
            if strcmpi(solname(k+4+ipos),'P'),break,end
        end
        ntor = str2double(solname(k+4:k+4+ipos-1));
        k = strfind(solname,'Ptor');
        ptor = str2double(solname(k+4));
        refElTor = createReferenceElement(0,(ptor+1)^2,[]);
        for ieq = fields2plot
            %figure(fig), fig = fig+1;
            %tiledlayout(row_subplot,col_subplot)
            clim1 = min(up(:,ieq));
            clim2 = max(up(:,ieq));
            for itor = 1:ntpos
                figure(fig), fig = fig+1;
                %subplot(row_subplot,col_subplot,itor)
                %nexttile
                upol = extractSolutionInAtGivenTheta(up(:,ieq),T,refEl,refElTor,tpos(itor));
                if diff
                    upol0 = extractSolutionInAtGivenTheta(up0(:,ieq),T,refEl,refElTor,tpos(itor));
                    upol0bis = extractSolutionInAtGivenTheta(up0(:,ieq),T,refEl,refElTor,tpos(itor));
                elseif average
                    upol0 = upm(:,ieq);
                    upol0bis = upm(:,ieq);
                else
                    upol0 = 0.;
                    upol0bis = 1.;
                end
                if ieq == 1
                    upol = (upol - upol0)*100./upol0bis;
                    %clim1 = -5; %min(upol);
                    %clim2 = 5; %max(upol);
                    title(['N, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'N';

                elseif ieq == 2
                    stop
                    upol = (upol - upol0);
                    clim1 = min(upol);
                    clim2 = max(upol);
                    if neq == 2
                        title(['Mach, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                        f = 'Mach';
                        
                    else
                        title(['u para, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                        f = 'u_para';
                    end
                elseif ieq == 3
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Ei, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Ei';
                elseif ieq == 4
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Ee, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Ee';
                elseif ieq == 5
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Pi, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Pi';
                elseif ieq == 6
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Pe, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Pe';
                elseif ieq == 7
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Ti, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Ti';
                elseif ieq == 8
                    upol = (upol - upol0)*100./upol0bis;
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Te, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Te';
                elseif ieq == 9
                    upol = (upol - upol0);
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Cs, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Cs';
                elseif ieq == 10
                    upol = (upol - upol0);
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['Mach, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'Mach';
                else
                    upol = (upol - upol0);
                    clim1 = min(upol);
                    clim2 = max(upol);
                    title(['???, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
                    f = 'bob';
                end
                plotSolution(X/lscale,T,upol,refEl);
                if mesh2plot
                    hold on, plotMesh(X/lscale,T,elemType,[]);
                end
                xlim([xlim1,xlim2]);
                ylim([ylim1,ylim2]);
                caxis([clim1,clim2]);
                box on;
                %caxissave(itor,:) = caxis;
                drawnow
                if saveplot
                    saveas(gcf, [savedir, datestr(now,'yyyy_mm'),'_',str_diff,'_',f,'_',pad(int2str(itor),3,'left','0'),'.png'])
                end
            end
        end
    end
elseif strcmpi(solname(4:5),'2D') %2D
    close all
    if consSol
        for ieq = fields2plot
            figure(fig), fig = fig+1;
            clim1 = min(u(:,ieq));
            clim2 = max(u(:,ieq));
            if ieq == 1
                title('N')
            elseif ieq == 2
                title('Gamma')
            elseif ieq == 3
                title('NEi')
            elseif ieq == 4
                title('NEe')
            else
                title('???')
            end
            plotSolution(X/lscale,T,u(:,ieq),refEl);
            if mesh2plot
                hold on, plotMesh(X/lscale,T,elemType,[]);
            end
            xlim([xlim1,xlim2]);
            ylim([ylim1,ylim2]);
            box on;
        end
    else
        close all
        for ieq = fields2plot
            clim1 = min(up(:,ieq));
            clim2 = max(up(:,ieq));
            figure(fig), fig = fig+1;
            if ieq == 1
                title('N')
            elseif ieq == 2
                if neq == 2
                    title('Mach')
                else
                    title('u para')
                end
            elseif ieq == 3
                title('Ei')
            elseif ieq == 4
                title('Ee')
            elseif ieq == 5
                title('Pi')
            elseif ieq == 6
                title('Pe')
            elseif ieq == 7
                title('Ti')
            elseif ieq == 8
                title('Te')
            elseif ieq == 9
                title('Cs')
            elseif ieq == 10
                title('Mach')
            else
                title('???')
            end
            plotSolution(X/lscale,T,up(:,ieq),refEl);
            if mesh2plot
                hold on, plotMesh(X/lscale,T,elemType,[]);
            end
            xlim([xlim1,xlim2]);
            ylim([ylim1,ylim2]);
            box on;
            %caxis([clim1,clim2]);
            %caxissave(itor,:) = caxis;
        end
    end
end

% %Plot solution
% for ieq = 1:neq
%     if strcmpi(solname(4:5),'2D') %2D
%         figure(fig), fig = fig+1;
%         if ieq~=2
%             if ieq==1
%                 plotSolution(X/lscale,T,u(:,ieq),refEl,5)
%             else
%                 plotSolution(X/lscale,T,u(:,ieq)./u(:,1),refEl,5)
%             end
%         else
%             plotSolution(X/lscale,T,u(:,ieq)./u(:,1),refEl,5)
%         end
%         if mesh2plot
%             hold on, plotMesh(X/lscale,T,elemType,[])
%         end
%         if ieq == 1
%             title('N')
%         elseif ieq == 2
%             title('u para')
%         elseif ieq == 3
%             title('Ei')
%         elseif ieq == 4
%             title('Ee')
%         else
%             title('???')
%         end
%     elseif strcmpi(solname(4:5),'3D') %3D
%         k = strfind(solname,'Ntor');
%         for ipos = 1:3
%             if strcmpi(solname(k+4+ipos),'P'),break,end
%         end
%         ntor = str2double(solname(k+4:k+4+ipos-1));
%         k = strfind(solname,'Ptor');
%         ptor = str2double(solname(k+4));
%         refElTor = createReferenceElement(0,(ptor+1)^2,[]);
%         %figure(fig), fig = fig+1;
%         %tiledlayout(row_subplot,col_subplot)
%         clim1 = min(u(:,ieq));
%         clim2 = max(u(:,ieq));
%         for itor = 1:ntpos
%             figure(fig), fig = fig+1;
%             %subplot(row_subplot,col_subplot,itor)
%             %nexttile
%             if ieq == 1
%                 title(['N, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
%             elseif ieq == 2
%                 title(['Gamma, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
%             elseif ieq == 3
%                 title(['NEi, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
%             elseif ieq == 4
%                 title(['NEe, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
%             else
%                 title(['???, plane :',num2str(tpos(itor)),'/',num2str(tmax)])
%             end
%             upol = extractSolutionInAtGivenTheta(u(:,ieq),T,refEl,refElTor,tpos(itor));
%             plotSolution(X/lscale,T,upol,refEl);
%             xlim([xlim1,xlim2]);
%             ylim([ylim1,ylim2]);
%             %caxis([clim1,clim2]);
%             %     caxissave(itor,:) = caxis;
%             drawnow
%         end
%     end
% end

function res = faceNodes_aux_quads(nOfElementNodes)
switch nOfElementNodes
    case 4 %P1
        res = [1 2; 2 3; 3 4; 4 1];
    case 9 %P2
        res = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
    case 16 %P3
        res = [1 5 6 2; 2 7 8 3; 3 9 10 4; 4 11 12 1];
    case 25 %P4
        res = [1 5 6  7 2; 2 8 9 10 3; 3 11 12 13 4; 4 14 15 16 1];
    case 36 %P5
        res = [1 5 6  7  8 2; 2 9 10 11 12 3; 3 13 14 15 16 4; 4 17 18 19 20 1];
    case 49 %P6
        res = [1 5 6  7  8 9 2; 2 10 11 12 13 14 3; 3 15 16 17 18 19  4; 4 20 21 22 23 24  1];
    case 64 %P7
        res = [1 5 6  7  8 9 10 2; 2 11 12 13 14 15 16 3; 3 17 18 19 20 21 22 4; 4 23 24 25 26 27 28 1];
    case 81 %P8
        res = [1 5 6  7  8 9 10 11 2; 2 12 13 14 15 16 17 18 3; 3 19 20 21 22 23 24 25 4; 4  26 27 28 29 30 31 32 1];
end
end

function res = faceNodes_aux(nOfElementNodes)
switch nOfElementNodes
    case 3 %P1
        res = [1 2; 2 3; 3 1];
    case 6 %P2
        res = [1 4 2; 2 5 3; 3 6 1];
    case 10 %P3
        res = [1 4 5 2; 2 6 7 3; 3 8 9 1];
    case 15 %P4
        res = [1 4 5 6 2; 2 7 8 9 3; 3 10 11 12 1];
    case 21 %P5
        res = [1 4:7 2; 2 8:11 3; 3 12:15 1];
    case 28 %P6
        res = [1 4:8 2; 2 9:13 3; 3 14:18 1];
    case 36 %P7
        res = [1 4:9 2; 2 10:15 3; 3 16:21 1];
    case 45 %P8
        res = [1 4:10 2; 2 11:17 3; 3 18:24 1];
    case 55 %P9
        res = [1 4:11 2; 2 12:19 3; 3 20:27 1];
    case 66 %P10
        res = [1 4:12 2; 2 13:21 3; 3 22:30 1];
    case 78 %P11
        res = [1 4:13 2; 2 14:23 3; 3 24:33 1];
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