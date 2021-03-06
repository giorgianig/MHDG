% plot Fortran solution
clear
close all

%**********************************
% Parallel/serial
%**********************************
parallel=1;
nproc=48;

%**********************************
% Plot options
%**********************************
plotPhys = 1; % Plot physical values
plotCons = 0; % Plot conservative variables
% Dimensional (1) or non-dimensional (0) plots
cons_dimensional_plots = 0; % conservative variables
phys_dimensional_plots = 0; % physical variables
nref = 1; % plot order

%**********************************
% 3D stuff
%**********************************
ntpos =2;
theta = 2*pi;        % Toroidal angle considered
dth = 1e-14;
tpos = linspace(0+dth,theta-dth,ntpos);

%**********************************
% Printing out
%**********************************
printout = 0; % Print out plots
path2save = '/home/giorgio/Dropbox/PostDoc_Marseille/Latex/NGammaVortPot/';


%**********************************
% Solution
%**********************************
solpath = '/home/giorgio/Dropbox/Fortran/Results/';
meshpath = '/home/giorgio/Dropbox/Fortran/MHDG_ref3.0/matlab/Meshes/Parallel/';
solname = 'Sol2D_CircLimAlign_Quads_Nel98784_P1_DPe0.500E+00';





%**********************************
% Working...
%**********************************

pos = strfind(solname,'_P');
for i=1:10
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end

% Make sure that I look only for 1 process if not parallel
if ~parallel
    nproc = 1;
end

for iproc = 1:nproc
    
    % Load results
    if parallel
        meshname = [solname(7:pos) '_' num2str(iproc) '_' num2str(nproc) '.h5'];
        solname_comp = [solname '_' num2str(iproc) '_' num2str(nproc) '.h5'];
    else
        meshname = [solname(7:pos), '.h5'];
        solname_comp = [solname, '.h5'];
    end
    HDF5load([solpath,solname_comp]);
    Mesh = HDF5load([meshpath,meshname]);
    
    % Number of equations
    Neq = simulation_parameters.Neq;
    
    % Conservative and physical variables
    uc = transpose(reshape(u,[Neq,numel(u)/Neq]));
    up = cons2phys(uc,simulation_parameters);
    
    % Reference element on plane
    if Mesh.elemType==1
        elemType=0;
    elseif Mesh.elemType==0
        elemType=1;
    end
    refEl = createReferenceElement(elemType,size(Mesh.T,2));
    
    iplot = 0;
    if strcmpi(solname(4:5),'2D')
        %**********************************
        % 2D plot...
        %**********************************
        if plotCons
            for i=1:size(uc,2)
                uplot = uc(:,i);
                if cons_dimensional_plots
                    uplot = uc(:,i)*simulation_parameters.adimensionalization.reference_values_conservative_variables(i);
                end
                iplot = iplot +1;
                figure(iplot),hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref);axis off                
                name = simulation_parameters.physics.conservative_variable_names{i};
                title(name)
                if printout
                    readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                end
            end
        end
        if plotPhys
            for i=1:size(up,2)
                uplot = up(:,i);
                if phys_dimensional_plots
                    uplot = up(:,i)*simulation_parameters.adimensionalization.reference_values_physical_variables(i);
                end
                iplot = iplot +1;
                figure(iplot),hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref);axis off
                if i==3,plotMesh(Mesh.X,Mesh.T,elemType),end
                name = simulation_parameters.physics.physical_variable_names{i};
                title(name)
                if printout
                    readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                end
            end
        end
    elseif strcmpi(solname(4:5),'3D')
        %**********************************
        % 3D plot...
        %**********************************
        k = strfind(solname,'Ntor');
        for ipos = 1:3
            if strcmpi(solname(k+4+ipos),'P'),break,end
        end
        ntor = str2double(solname(k+4:k+4+ipos-1));
        k = strfind(solname,'Ptor');
        ptor = str2double(solname(k+4));
        refElTor = createReferenceElement(0,(ptor+1)^2,[]);
        for itor = 1:ntpos
            if plotCons
                for i=1:size(uc,2)
                    uplot = uc(:,i);
                    if cons_dimensional_plots
                        uplot = uc(:,i)*simulation_parameters.adimensionalization.reference_values_conservative_variables(i);
                    end
                    uplot = extractSolutionInAtGivenTheta(uplot,T,refEl,refElTor,tpos(itor));
                    iplot = iplot +1;
                    figure(iplot),hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref);axis off
                    name = simulation_parameters.physics.conservative_variable_names{i};
                    title([name,'- Plane ' num2str(itor)])
                    if printout
                        readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                    end
                    drawnow
                end
            end
            if plotPhys
                for i=1:size(up,2)
                    uplot = up(:,i);
                    if phys_dimensional_plots
                        uplot = up(:,i)*simulation_parameters.adimensionalization.reference_values_physical_variables(i);
                    end
                    uplot = extractSolutionInAtGivenTheta(uplot,T,refEl,refElTor,tpos(itor));
                    iplot = iplot +1;
                    figure(iplot),hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref);axis off
                    name = simulation_parameters.physics.physical_variable_names{i};
                    title([name,'- Plane ' num2str(itor)])
                    if printout
                        readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                    end
                    drawnow
                end
            end
        end
    end
end
% Mesh.maxx = max((X(:,1)+3.4)/lscale);
% Mesh.minx = min((X(:,1)+3.4)/lscale);
%
% testcase.n=60;
% testcase.xc=0;
% testcase.yc=0;
% b = defineMagneticField([X(T',1)+3.4,X(T',2)]/lscale);
% qres = transpose(reshape(q,[neq*2,numel(q)/2/neq]));
%
% Gradpar_phi = qres(:,7).*b(:,1)+qres(:,8).*b(:,2);
% Gradpar_n = qres(:,1).*b(:,1)+qres(:,2).*b(:,2);
% figure, plotSolution(X/lscale,T,Gradpar_phi-Mref*Gradpar_n./u(1:neq:end-3),refEl,5);axis off,title('Parallel courrent')
% if printout
%     readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,['ParallelCurrent',num2str(testnumber)])
% end
% Gradper_phix = qres(:,7)-Gradpar_phi.*b(:,1);
% Gradper_phiy = qres(:,8)-Gradpar_phi.*b(:,2);
%
% Gradper_nx = qres(:,1)-Gradpar_n.*b(:,1);
% Gradper_ny = qres(:,2)-Gradpar_n.*b(:,2);
%
% polacourr = sqrt( (Gradper_nx+Mref*Gradper_phix./u(1:neq:end-3)).^2+(Gradper_ny+Mref*Gradper_phiy./u(1:neq:end-3)).^2);
% figure, plotSolution(X/lscale,T,polacourr,refEl,5);axis off,title('Pola courrent')
% if printout
%     readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,['Polacourr',num2str(testnumber)])
% end
% figure, plotSolution(X/lscale,T,Gradper_phiy-Mref*Gradper_phiy./u(1:neq:end-3),refEl,5);axis off,title('Pola courrent y')
% readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,'PerpCurrenty')




