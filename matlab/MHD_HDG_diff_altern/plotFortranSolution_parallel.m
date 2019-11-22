% plot Fortran solution parallel
clear
close all

global Mesh theta ntor
kT=25;
lscale =  1.901*1e-3;
nproc = 4;

% for 3D
ntpos =2;
theta = 6.28;        % Toroidal angle considered

dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);


ntpos=1;
tpos = tpos(1);

%% West
% matmeshpath = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/';
% solpath = '/home/giorgio/Saves_MHDG_Marconi/West/ApplyingThreshold/Deuterium/';
% meshpath = '/home/giorgio/Saves_MHDG_Marconi/West/';
% solution = 'Sol_West_Hugo_h0.06_refCorn0.008_refCornCrit0.002_refSep0.01_P6_Diff.80000E-02.h5';
% solution = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_P8_Diff.86067E-01';
% solution = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_P8_Diff.38000E-01';
% solution = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_P8_Diff.16778E-01';
% solution = 'Sol_West_Hugo_h0.02_refCorn0.001_refSep0.01_P8_Diff.74076E-02_0001';

%% Goldston
matmeshpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
% solpath = '/home/giorgio/Goldston/provv/';

solution = 'Sol3D_Circle_lim_h02_P4_Ntor2Ptor1_DPe0.380E+00_NR0000';
% solution = 'Sol2D_Circle_lim_h02_P4_DPe0.380E+00_NR0001';

meshpath =  matmeshpath;


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
    
    
    
    if strcmpi(solution(4:5),'2D')
        figure()
        for ip = 1:nproc
            meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
            solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
            HDF5load([meshpath,meshname])
            HDF5load([solpath,solname])
            u = u';
            if elemType==0
                eltype=1;
            else 
                eltype=0;
            end
            refEl = createReferenceElement(1,size(T,2));
            plotSolution(X/lscale,T,u(1:2:end-1),refEl,5)
            hold on, plotMesh(X/lscale,T,eltype,[])
        end
    elseif strcmpi(solution(4:5),'3D')
        k = strfind(solution,'Ntor');
        for ipos = 1:3
            if strcmpi(solution(k+4+ipos),'P'),break,end
        end
        ntor = str2double(solution(k+4:k+4+ipos-1));
        k = strfind(solution,'Ptor');
        ptor = str2double(solution(k+4));
        refElTor = createReferenceElement(0,(ptor+1)^2,[]);
        
        for itor = 1:ntpos
            figure(),clf
            for ip = 1:nproc
                meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
                solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
                HDF5load([meshpath,meshname])
                HDF5load([solpath,solname])
                u = u';
                refEl = createReferenceElement(1,size(T,2));
                upol = extractSolutionInAtGivenTheta(u(1:2:end-1),T,refEl,refElTor,tpos(itor));
                
                plotSolution(X/lscale,T,upol,refEl)
                %     caxissave(itor,:) = caxis;
                %         title(['Plane ' num2str(itor)])
                drawnow
            end
        end
        
        
    end
    title('Density')
    %
    %     figure
    %     for ip = 1:nproc
    %         meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    %         solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
    %         HDF5load([meshpath,meshname])
    %         HDF5load([solpath,solname])
    %         u = u';
    %         refEl = createReferenceElement(1,size(T,2));
    %         plotSolution(X/lscale,T,u(2:2:end),refEl,5);
    %         hold on
    %     end
    %     title('Gamma')
    
%     figure
%     for ip = 1:nproc
%         meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
%         solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
%         HDF5load([meshpath,meshname])
%         HDF5load([solpath,solname])
%         
%         %        X = X/0.75*0.287;
%         
%         
%         u = u';
%         refEl = createReferenceElement(1,size(T,2));
%         plotSolution(X/lscale,T,u(2:2:end)./u(1:2:end-1)./sqrt(kT),refEl,5);
%         axis off
%         hold on
%     end
%     
%     %     colorbar off
%     %     caxis([ -2.    2.5])
%     %     axis([ -0.6718   -0.5587    0.3895    0.4787])
%     title('Mach')
    %
    %     figure
    %     for ip = 1:nproc
    %         meshname = [solution(5:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    %         solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
    %         HDF5load([meshpath,meshname])
    %         HDF5load([solpath,solname])
    %         u = col(scdiff_nodes');
    %
    %         r = sqrt(X(T',1).^2+X(T',2).^2);
    % %         u(r<0.25) = 0;
    %
    %         refEl = createReferenceElement(1,size(T,2));
    %         plotSolution(X/lscale,T,u,refEl,5);
    %         hold on
    %     end
    %     axis off
    % %     axis([ -0.6718   -0.5587    0.3895    0.4787])
    %     title('Shock capturing coefficient')
    
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