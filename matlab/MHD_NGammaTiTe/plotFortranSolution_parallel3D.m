% plot Fortran solution parallel
clear
close all

global Mesh Mref neq ntor theta
Mref = 12.5;
lscale =  1.901*1e-3;
nproc = 8;
neq = 4;
nref = 5;

elemType = 1; %1: Triangle 0:Quads

ntpos = 17; % number of toroidal positions to interpolate from the simulation
theta = 2*pi/9; % Toroidal angle in the simulation (I don't know why it is called theta but do not change the name)
dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);
cons_phys = true; % true to plot physical values (N,u_para, T,etc)
field2plot = 10; % if cons_phys  true [1:N, 2:u_para, 3:Ei, 4:Ee, 5:Pi, 6:Pe, 7:Ti, 8:Te, 9:Cs, 10:M] else [1:N, 2:Gamma_para, 3:NEi, 4:NEe]
average = true; % true if toroidal average


matmeshpath = '/home/bluce/MHDG_sim/West10861_P4/3D/NGT/parall/00/ref/ripple/';
solpath = '/home/bluce/MHDG_sim/West10861_P4/3D/NGT/parall/00/ref/ripple/ntor1ptor4/';
solution = 'Sol3D_West_YesHole_Nel10861_P4_Ntor1Ptor4_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';
meshpath =  '/home/bluce/MHDG_sim/West10861_P4/3D/NGT/parall/00/ref/ripple/';


%% start
cases =1;
if cons_phys
    legend = ["N" "u_para" "Ei" "Ee" "Pi" "Pe" "Ti" "Te" "Cs" "M"];
else
    legend = ["N" "Gamma" "NEi" "NEe"];
end


Mesh.lscale = lscale;
pos = strfind(solution,'_P');
for i=1:10
    if strcmp(solution(pos+i),'_')
        pos = pos+i-1;
        break
    end
end

if average
    figure
else
    for itpos = 1:ntpos
        figure(itpos)
    end
end

for ip = 1:nproc
    meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    solname = [solution '_ip' num2str(ip) '_it1_np' num2str(nproc) '.h5'];
    HDF5load([meshpath,meshname])
    HDF5load([solpath,solname])
    u = transpose(reshape(u,[neq,numel(u)/neq]));
    if cons_phys
        up = cons2phys(u);
    end
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
    if average
        sz = size(extractSolutionInAtGivenTheta(u(:,1),T,refEl,refElTor,tpos(1)),1);
        um = zeros(sz,ntpos);
        for itpos = 1:ntpos
            if cons_phys
                um(:,itpos) = extractSolutionInAtGivenTheta(up(:,field2plot),T,refEl,refElTor,tpos(itpos));
            else
                um(:,itpos) = extractSolutionInAtGivenTheta(u(:,field2plot),T,refEl,refElTor,tpos(itpos));
            end
        end
        u2D = mean(um,2);
        plotSolution(X/lscale,T,u2D,refEl,nref);
        hold on
        axis off,title(strcat('<', legend(field2plot), '> (\phi)'))
    else
        for itpos = 1:ntpos
            if cons_phys
                u2D = extractSolutionInAtGivenTheta(up(:,field2plot),T,refEl,refElTor,tpos(itpos));
            else
                u2D = extractSolutionInAtGivenTheta(u(:,field2plot),T,refEl,refElTor,tpos(itpos));
            end
            figure(itpos)
            plotSolution(X/lscale,T,u2D,refEl,nref);
            hold on
            [num,dem]=rat(tpos(itpos)/pi);
            if dem >1e5               
                axis off,title(strcat(legend(field2plot),' \phi=0'))
            else
                axis off,title(strcat(legend(field2plot),' \phi=',num2str(num),'\pi','/',num2str(dem)))
            end
        end
    end
end


% load([ matmeshpath solution(5:pos) '.mat' ])
% plotMesh(X,T)
% clear Tb
% aux = whos('Tb*');
% boundaryNames = {aux.name};
% clear aux
% nOfBound = numel(boundaryNames);
% for ibound = 1:nOfBound
%     name = boundaryNames{ibound};
%     eval(['aux_bound = Tb_' name(4:end) ';'])
%     for ifa = 1:size(aux_bound,1)
%         face = aux_bound(ifa,:);
%         plot(X(face,1),X(face,2),'k-','linewidth',1)
%         hold on
%     end
% end
% readyforprintjpeg([8 6],24,[],[],[],[],[],'GoldMach1e-3_NoDrift_caxis-auto')