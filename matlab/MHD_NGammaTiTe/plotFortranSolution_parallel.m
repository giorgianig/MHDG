% plot Fortran solution parallel
clear
close all

global Mesh
lscale =  1.901*1e-3;
nproc = 4;
neq = 4;
nref = 5;


matmeshpath = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';
solution = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P4_DPe0.380E+00_DPai0.300E+06_DPae0.100E+08_0001';
meshpath =  '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';


%% start
cases =1;

for icase = 1:numel(cases)
    
    
    solution = [solution(1:end-1) num2str(cases(icase))];
    
    Mesh.lscale = lscale;
    pos = strfind(solution,'_P');
    for i=1:10
        if strcmp(solution(pos+i),'_')
            pos = pos+i-1;
            break
        end
    end
  
    figure
    for ip = 1:nproc
        meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
        solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
        HDF5load([meshpath,meshname])
        HDF5load([solpath,solname])
        u = transpose(reshape(u,[neq,numel(u)/neq]));
        refEl = createReferenceElement(1,size(T,2));
            plotSolution(X/lscale,T,u(:,1),refEl,nref);
        hold on
    end
axis off,title('U1')

    figure
    for ip = 1:nproc
        meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
        solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
        HDF5load([meshpath,meshname])
        HDF5load([solpath,solname])
        u = transpose(reshape(u,[neq,numel(u)/neq]));
        refEl = createReferenceElement(1,size(T,2));
            plotSolution(X/lscale,T,u(:,2),refEl,nref);
        hold on
    end
axis off,title('U2')

    figure
    for ip = 1:nproc
        meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
        solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
        HDF5load([meshpath,meshname])
        HDF5load([solpath,solname])
        u = transpose(reshape(u,[neq,numel(u)/neq]));
        refEl = createReferenceElement(1,size(T,2));
            plotSolution(X/lscale,T,u(:,3),refEl,nref);
        hold on
    end
axis off,title('U3')

    figure
    for ip = 1:nproc
        meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
        solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
        HDF5load([meshpath,meshname])
        HDF5load([solpath,solname])
        u = transpose(reshape(u,[neq,numel(u)/neq]));
        refEl = createReferenceElement(1,size(T,2));
            plotSolution(X/lscale,T,u(:,4),refEl,nref);
        hold on
    end
axis off,title('U4') 
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