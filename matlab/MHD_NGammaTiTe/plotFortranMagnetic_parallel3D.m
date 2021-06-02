% plot Fortran solution parallel
clear
close all

global Mesh Mref neq ntor theta
Mref = 12.5;
lscale =  1.901*1e-3;
nproc = 4;
neq = 4;
nref = 5;

R0 = 3.4; %2.5; %3.4;

elemType = 1; %1: Triangle 0:Quads for matlab function but be careful, Fortran use the opposite!!!!

ntpos = 17; % number of toroidal positions to interpolate from the simulation
theta = 2*pi/9; % Toroidal angle in the simulation (I don't know why it is called theta but do not change the name)
dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);
average = false; % true if toroidal average

rmp=false;
ripple=true;


meshpath = '/home/bluce/MHDG_sim/InfLim909_P5/3D/NG/parall/00/ripple/';
solpath = '/home/bluce/MHDG_sim/InfLim909_P5/3D/NG/parall/00/ripple/';
solution = 'Sol3D_Circ_InfLIM_Quads_YesHole_Nel909_P5_Ntor3Ptor4_DPe0.100E+01';


%% start
cases =1;

Mesh.lscale = lscale;
pos = strfind(solution,'_P');
for i=1:10
    if strcmp(solution(pos+i),'_')
        pos = pos+i-1;
        break
    end
end

figure % Plot of the coils in 3D
for ip = 1:nproc
    meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    solname = [solution '_ip' num2str(ip) '_it1_np' num2str(nproc) '.h5'];
    HDF5load([meshpath,meshname])
    HDF5load([solpath,solname])
    
    % Correcting mess up with refElem
    if elemType==1
        elemType=0; %quads
    elseif elemType==0
        elemType=1; %triangle
    end
    
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
    nt = size(magnetic_field,1)/Nnodes; % toroidal number of mesh to plot
    smat = size(X,1)*nt;
    X3D = zeros(smat,3);
    T3D = zeros(size(T,1)*nt,size(T,2));
    
    r = sqrt((X(:,1) - R0).^2 + X(:,2).^2);
    Z = X(:,2);
    R = r.*cos(asin(Z./r));
    
    for i=1:nt
        tor_angle = theta*i/nt;
        X3D((i-1)*size(X,1)+1:i*size(X,1),1) = R.*cos(tor_angle);
        X3D((i-1)*size(X,1)+1:i*size(X,1),2) = R.*sin(tor_angle);
        X3D((i-1)*size(X,1)+1:i*size(X,1),3) = Z;
        T3D((i-1)*size(T,1)+1:i*size(T,1),:) = T + (i-1)*size(X,1);
    end
    %Connectivity for the faces
    patchFaces = T3D(:,oFaceNodes);
    %Plot 3D mesh
    patchHandle = patch('Faces',patchFaces,'Vertices',X3D,'FaceColor','none','EdgeAlpha',1,'linewidth',0.4, 'linestyle',':');
    
    % Plot of the coils
    % Load coils
    if ripple
        n = numel(coils_ripple(:,1)); % 2 positions per line (3 coordinates)
        coils_ripple = coils_ripple.*lscale;
        for i=1:n % segments
            line([coils_ripple(i,1) coils_ripple(i,2)],[coils_ripple(i,3) coils_ripple(i,4)],[coils_ripple(i,5) coils_ripple(i,6)], 'color', 'blue', 'linewidth',2);
        end
    end
    if rmp
        tmp = size(coils_rmp);
        nbRow = tmp(1);
        coils_rmp = coils_rmp.*lscale;
        for l=1:nbRow
            n = numel(coils_rmp(1,1,:)); % 2 positions per line (3 coordinates)
            for i=1:n % segments
                line([coils_rmp(l,1,i) coils_rmp(l,2,i)],[coils_rmp(l,3,i) coils_rmp(l,4,i)],[coils_rmp(l,5,i) coils_rmp(l,6,i)], 'color', 'red', 'linewidth',2);
            end
        end
    end

    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    box on;
    hold on;
end
% stop
% Plot magnetic fields (perturbation)
dimB = size(magnetic_field,1);
ntposB = dimB/Nnodes;
tposB = linspace(0+dth,theta-dth,ntposB);
for itpos = 1:ntposB
    figure(itpos)
end

for ip = 1:nproc
    meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    solname = [solution '_ip' num2str(ip) '_it1_np' num2str(nproc) '.h5'];
    HDF5load([meshpath,meshname])
    HDF5load([solpath,solname])   
    
    BR = magnetic_field(:,1);
    BZ = magnetic_field(:,2);
    BPhi = magnetic_field(:,3);
    Bmod = sqrt(BR.^2 + BZ.^2 + BPhi.^2);

    BRpert = magnetic_perturbation(:,1);
    BZpert = magnetic_perturbation(:,2);
    BPhipert = magnetic_perturbation(:,3);
    Bmodpert = sqrt(BRpert.^2 + BZpert.^2 + BPhipert.^2);

    % Correcting mess up with refElem
    if elemType==1
        elemType=0;
    elseif elemType==0
        elemType=1;
    end
    % Creating reference element (2D)
    refEl = createReferenceElement(elemType,size(T,2));
    % Creating reference element (3D)
    k = strfind(solname,'Ntor');
    for ipos = 1:3
        if strcmpi(solname(k+4+ipos),'P'),break,end
    end
    ntor = str2double(solname(k+4:k+4+ipos-1));
    k = strfind(solname,'Ptor');
    ptor = str2double(solname(k+4));
    refElTor = createReferenceElement(0,(ptor+1)^2,[]);
    
    dimB = size(magnetic_field,1);
    ntposB = dimB/Nnodes;
    tposB = linspace(0+dth,theta-dth,ntposB);
    clim1 = 0;
    clim2 = 1e-3;
    
    if average
        BRmean = zeros(Nnodes,ntposB);
        BZmean = zeros(Nnodes,ntposB);
        BPhimean = zeros(Nnodes,ntposB);
        for itpos = 1:ntposB
            BRmean(:,itpos) = BR((itpos-1)*Nnodes+1:itpos*Nnodes);
            BZmean(:,itpos) = BZ((itpos-1)*Nnodes+1:itpos*Nnodes);
            BPhimean(:,itpos) = BPhi((itpos-1)*Nnodes+1:itpos*Nnodes);
        end
        BRmean = mean(BRmean,2);
        BZmean = mean(BZmean,2);
        BPhimean = mean(BPhimean,2);
        for itpos = 1:ntposB
            sol = sqrt((BR((itpos-1)*Nnodes+1:itpos*Nnodes)-BRmean).^2 + (BZ((itpos-1)*Nnodes+1:itpos*Nnodes)-BZmean).^2 + (BPhi((itpos-1)*Nnodes+1:itpos*Nnodes)-BPhimean).^2);
%             sol = sqrt((BRmean).^2 + (BZmean).^2 + (BPhimean).^2);
            figure(itpos)
            plotSolution(X/lscale,T,sol,refEl,nref);
            hold on
            [num,dem]=rat(tposB(itpos)/pi);
            if dem >1e5
                axis off,title(strcat('|Bpert|av, plane :',' \phi=0'))
            else
                axis off,title(strcat('|Bpert|av, plane :',' \phi=',num2str(num),'\pi','/',num2str(dem)))
            end
            caxis([clim1, clim2]);
            box on;
        end
    else
        for itpos = 1:ntposB
            figure(itpos)
            plotSolution(X/lscale,T,Bmodpert((itpos-1)*Nnodes+1:itpos*Nnodes),refEl,nref);
            hold on
            [num,dem]=rat(tposB(itpos)/pi);
            if dem >1e5
                axis off,title(strcat('|Bpert|, plane :',' \phi=0'))
            else
                axis off,title(strcat('|Bpert|, plane :',' \phi=',num2str(num),'\pi','/',num2str(dem)))
            end
            caxis([clim1, clim2]);
            box on;
        end
    end
end




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