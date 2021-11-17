% plot Fortran solution
clear
close all

global ntor theta

%**********************************
% Parallel/serial
%**********************************
parallel=0; % true=1
nproc=8; % number of MPI task for parallel (mesh separation)

%**********************************
% Plot options
%**********************************
plotPhys = 1; % Plot physical values (true=1)
plotCons = 0; % Plot conservative variables
% Dimensional (1) or non-dimensional (0) plots
phys_dimensional_plots = 0; % physical variables
cons_dimensional_plots = 0; % conservative variables
nref = 3; % plot order
startPlot = 1; %default: 1 Starting number for plot (useful if no close all)
gatherplot = 0; %True to gather everything in one figure
cbound = 0; %True: bound colorbar axis (to adapt the boundaries see below)

%**********************************
% 3D stuff
%**********************************
ntpos =3; % Number of plot along the toroidal direction
dth = 1e-14; % starting point
average = 0; % true if toroidal average (for magnetic field (ripple,rmp), plot |B -<B>| instead of |Bpert|)
rmp = 0; % true if rmp
ripple = 0; % true if ripple

%**********************************
% Printing out
%**********************************
printout = 0; % Print out plots
path2save = 'Img_MHDG/';

%**********************************
% Solution
%**********************************
HOME = '/home/giorgio/Dropbox/Fortran/MHDG_ref3.0/test/';

% solpath = [HOME, 'West10861_P4/2D/NGT/parall/00/div16/2021_06_08_old_Bohmbc_02/'];
% meshpath = [HOME, 'West10861_P4/2D/NGT/parall/00/div16/'];
% solname = 'Sol2D_West_YesHole_Nel10861_P4_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'West10861_P4/2D/NGT/parall/00/div16/2021_07_06_Neutral/'];
% meshpath = [HOME, 'West10861_P4/2D/NGT/parall/00/div16/'];
% solname = 'Sol2D_West_YesHole_Nel10861_P4_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath =[HOME, 'test/NGT3D_div8/'];
% meshpath = [HOME, 'test/NGT3D_div8/'];
% solname = 'Sol3D_CircLimAlign_Quads_Nel208_P4_Ntor1Ptor5_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'West3541_P4/2D/NGT/parall/00/div16/2021_06_21_NGT/'];
% meshpath = [HOME, 'West3541_P4/2D/NGT/parall/00/div16/'];
% solname = 'Sol2D_West_YesHole_Nel3541_P4_DPe0.200E+02_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'West3541_P4/3D/NGT/parall/00/ref/ripple/div16/'];
% meshpath = [HOME, 'West3541_P4/3D/NGT/parall/00/ref/ripple/div16/'];
% solname = 'Sol3D_West_YesHole_Nel3541_P4_Ntor1Ptor4_DPe0.138E+02_DPai0.314E+06_DPae0.105E+08_NR0006';

% solpath = [HOME, 'InfLim909_P5/2D/NGT/parall/00/'];
% meshpath = [HOME, 'InfLim909_P5/2D/NGT/parall/00/';
% solname = 'Sol2D_Circ_InfLIM_Quads_YesHole_Nel909_P5_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'test/NGT3D_div8/'];
% meshpath = [HOME, 'test/NGT3D_div8/'];
% solname = 'Sol3D_CircLimAlign_Quads_Nel208_P4_Ntor1Ptor5_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'test/NGTN2D_div8/'];
% meshpath = [HOME, 'test/NGTN2D_div8/'];
% solname = 'Sol2D_CircLimAlign_Quads_Nel208_P4_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

% solpath = [HOME, 'West1902_P6/2021_07_06_NGTN2D/'];
% meshpath = [HOME, 'West1902_P6/2021_07_06_NGTN2D/'];
% solname = 'Sol2D_West_YesHole_Nel1902_P6_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08';

solpath = [HOME];
meshpath = [HOME];
solname = 'Sol2D_ITER_YesHoleSmooth_Quads_Nel17611_P4_DPe0.308E+01_DPai0.314E+06_DPae0.105E+08_NR0005';



%**********************************
% Bound axis
%**********************************
if plotPhys && cbound
    if phys_dimensional_plots
        % rho, u, Ei, Ee, Pi, Pe, Ti, Te, Csi, M , rhon
        cpaxis= [[0 1.2e19]; [-0.25 0.25]; [0 3.5e9]; [0 3.5e9]; [0 2.3e28]; [0 2.3e28]; [0 50]; [0 50]; [2 7e4]; [-1.1 1.1]; [0 4.0e-3]];
    else
        % rho, u, Ei, Ee, Pi, Pe, Ti, Te, Csi, M , rhon
        cpaxis= [[0 1.2]; [-0.25 0.25]; [0 20]; [0 20]; [0 1]; [0 1]; [0 1]; [0 1]; [1 5]; [-1.1 1.1]; [0 4.0e-3]];
    end
else
    cpaxis= [[0 1.2]; [-0.25 0.25]; [0 20]; [0 20]; [0 1]; [0 1]; [0 1]; [0 1]; [1 5]; [-1.1 1.1]; [0 4.0e-3]];
end
if plotCons && cbound
    if cons_dimensional_plots
        % rho, Gamma (Nu), NEi, NEe, rhon
        ccaxis= [[0 1.2e19]; [-1.6e16 1.6e16]; [0 3.5e28]; [0 3.5e28]; [0 3.0e-3]];
    else
        ccaxis= [[0 1.2]; [-0.9 0.9]; [0 20]; [0 20]; [0 3.0e-3]];
    end
else
    ccaxis= [[0 1.2]; [-0.9 0.9]; [0 20]; [0 20]; [0 3.0e-3]];
end

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
        if strcmpi(solname(4:5),'3D')
            solname_comp = [solname '_ip' num2str(iproc) '_it1_np' num2str(nproc) '.h5'];
        else
            solname_comp = [solname '_' num2str(iproc) '_' num2str(nproc) '.h5'];
        end
    else
        meshname = [solname(7:pos), '.h5'];
        solname_comp = [solname, '.h5'];
    end
    HDF5load([solpath,solname_comp]);
    Mesh = HDF5load([meshpath,meshname]);
    Mesh.lscale = simulation_parameters.adimensionalization.length_scale;
    % Number of equations
    Neq = simulation_parameters.Neq;
    
    % Conservative and physical variables
    uc = transpose(reshape(u,[Neq,numel(u)/Neq]));
    up = cons2phys(uc,simulation_parameters);
    nc = size(uc,2)+1;
    np = size(up,2)+1;
    [nrowc, ncolc, nrowp, ncolp] = setSubplot(simulation_parameters);
    
    % Reference element on plane
    if Mesh.elemType==1
        elemType=0;
    elseif Mesh.elemType==0
        elemType=1;
    end
    refEl = createReferenceElement(elemType,size(Mesh.T,2));
    
    if cbound
        cont = 2;
    else
        cont = 0;
    end
    
    iplot = startPlot;
    if strcmpi(solname(4:5),'2D')
        %**********************************
        % 2D plot...
        %**********************************
        if iproc == 1 && gatherplot
            if plotCons
                fc = figure('Name','Conservative values','WindowState','maximized');
                hac = tight_subplot(nrowc,ncolc,[.01 .03],[.2 .03],[.03 .03]);
            end
            if plotPhys
                fp = figure('Name','Physical values','WindowState','maximized');
                hap = tight_subplot(nrowp,ncolp,[.01 .03],[.2 .03],[.03 .03]);
            end
        end
        if gatherplot
            if plotCons
                icplot = 0; icplot = icplot+1; figure(fc), axes(hac(icplot));
                hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                name = 'Mesh'; title(name), axis equal
            end
            if plotPhys
                ipplot = 0; ipplot = ipplot+1; figure(fp), axes(hap(ipplot));
                hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                name = 'Mesh'; title(name), axis equal
            end
        else
            figure(iplot)
            hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
            name = 'Mesh'; title(name)
            if printout
                readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
            end
        end
        if plotCons
            for ii=1:size(uc,2)
                uplot = uc(:,ii);
                if cons_dimensional_plots
                    uplot = uc(:,ii)*simulation_parameters.adimensionalization.reference_values_conservative_variables(ii);
                end
                if gatherplot
                    icplot = icplot +1; figure(fc), axes(hac(icplot));
                else
                    iplot = iplot +1; figure(iplot)
                end
                name = simulation_parameters.physics.conservative_variable_names{ii};
                if strcmpi(name,'rhon')
                    hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,uplot,refEl,nref,cont,1,0,ccaxis(ii,:));axis off
                else
                    hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,uplot,refEl,nref,cont,0,0,ccaxis(ii,:));axis off
                end
                
                title(name), axis equal
                if printout && ~gatherplot
                    readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                end
            end
        end
        if plotPhys
            for ii=1:size(up,2)
                uplot = up(:,ii);
                if phys_dimensional_plots
                    uplot = up(:,ii)*simulation_parameters.adimensionalization.reference_values_physical_variables(ii);
                end
                if gatherplot
                    ipplot = ipplot +1; figure(fp), axes(hap(ipplot)); 
                else
                    iplot = iplot +1; figure(iplot)
                end
                name = simulation_parameters.physics.physical_variable_names{ii};
                if strcmpi(name(1:4),'rhon')
                hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,abs(uplot)+1e-9,refEl,nref,cont,1,0,cpaxis(ii,:));axis off
                caxis([-4 0])
                else
                hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,uplot,refEl,nref,cont,0,0,cpaxis(ii,:));axis off
                end
                name = simulation_parameters.physics.physical_variable_names{ii};
                title(name), axis equal
                if printout && ~gatherplot
                    readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                end
            end
        end
    elseif strcmpi(solname(4:5),'3D')
        %**********************************
        % 3D plot...
        %**********************************
        %Toroidal angle considered
        theta = simulation_parameters.numerics.Max_extention_in_the_toroidal_direction;
        tpos = linspace(0+dth,theta-dth,ntpos);
        k = strfind(solname,'Ntor');
        for ipos = 1:3
            if strcmpi(solname(k+4+ipos),'P'),break,end
        end
        ntor = str2double(solname(k+4:k+4+ipos-1));
        k = strfind(solname,'Ptor');
        ptor = str2double(solname(k+4));
        refElTor = createReferenceElement(0,(ptor+1)^2,[]);
        if average
            if iproc == 1 && gatherplot
                if plotCons
                    fc = figure('Name','Conservative values','WindowState','maximized');
                    hac = tight_subplot(nrowc,ncolc,[.01 .03],[.2 .03],[.03 .03]);
                end
                if plotPhys
                    fp = figure('Name','Physical values','WindowState','maximized');
                    hap = tight_subplot(nrowp,ncolp,[.01 .03],[.2 .03],[.03 .03]);
                end
            end
            if gatherplot
                if plotCons
                    icplot = 0;icplot = icplot +1;figure(fc), axes(hac(icplot));
                    hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                    name = 'Mesh';
                    title(name), axis equal
                end
                if plotPhys
                    ipplot = 0; ipplot = ipplot +1;figure(fp), axes(hap(ipplot));
                    hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                    name = 'Mesh';
                    title(name), axis equal
                end
            else
                figure(iplot)
                hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                name = 'Mesh'; title(name)
                if printout
                    readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                end
            end
            if plotCons
                szc = size(extractSolutionInAtGivenTheta(uc(:,1),Mesh.T,refEl,refElTor,tpos(1)),1);
                umc = zeros(szc,ntpos,size(uc,2));
            end
            if plotPhys
                szp = size(extractSolutionInAtGivenTheta(up(:,1),Mesh.T,refEl,refElTor,tpos(1)),1);
                ump = zeros(szp,ntpos,size(up,2));
            end
            for itor = 1:ntpos
                if plotCons
                    for ii=1:size(uc,2)
                        umc(:,itor,ii) = extractSolutionInAtGivenTheta(uc(:,ii),Mesh.T,refEl,refElTor,tpos(itor));
                        if cons_dimensional_plots
                            umc(:,itor,ii) = umc(:,itor,ii)*simulation_parameters.adimensionalization.reference_values_conservative_variables(ii);
                        end
                    end
                end
                if plotPhys
                    for ii=1:size(up,2)
                        ump(:,itor,ii) = extractSolutionInAtGivenTheta(up(:,ii),Mesh.T,refEl,refElTor,tpos(itor));
                        if phys_dimensional_plots
                            ump(:,itor,ii) = ump(:,itor,ii)*simulation_parameters.adimensionalization.reference_values_physical_variables(ii);
                        end
                    end
                end
            end
            if plotCons
                for ii=1:size(uc,2)
                    uplot = mean(umc(:,:,ii),2);
                    if gatherplot
                        icplot = icplot +1; figure(fc), axes(hac(icplot)); 
                    else
                        iplot = iplot +1; figure(iplot)
                    end
                    hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,uplot,refEl,nref,cont,0,0,ccaxis(ii,:));axis off
                    name = simulation_parameters.physics.conservative_variable_names{ii};
                    title(strcat('<',name, '> (\phi)')), axis equal
                    if printout && ~gatherplot
                        readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                    end
                    drawnow
                end
                umc = 0;
            end
            if plotPhys
                for ii=1:size(up,2)
                    uplot = mean(ump(:,:,ii),2);
                    if gatherplot
                        ipplot = ipplot +1; figure(fp), axes(hap(ipplot));
                    else
                        iplot = iplot +1; figure(iplot)
                    end
                    hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,uplot,refEl,nref,cont,0,0,cpaxis(ii,:));axis off
                    name = simulation_parameters.physics.physical_variable_names{ii};
                    title(strcat('<',name, '> (\phi)')), axis equal
                    if printout
                        readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                    end
                    drawnow
                end
                ump = 0;
            end
        else
            for itor = 1:ntpos
                if iproc == 1 && gatherplot
                    if plotCons
                        fc(itor) = figure('Name','Conservative values','WindowState','maximized');
                        hac(:,itor) = tight_subplot(nrowc,ncolc,[.01 .03],[.2 .03],[.03 .03]);
                    end
                    if plotPhys
                        fp(itor) = figure('Name','Physical values','WindowState','maximized');
                        hap(:,itor) = tight_subplot(nrowp,ncolp,[.01 .03],[.2 .03],[.03 .03]);
                        ipplot = 0;
                    end
                end
                if gatherplot
                    if plotCons
                        icplot = 0;icplot = icplot +1;figure(fc(itor)), axes(hac(icplot,itor));
                        hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                        name = 'Mesh';
                        title(name), axis equal
                    end
                    if plotPhys
                        ipplot = 0; ipplot = ipplot +1;figure(fp(itor)), axes(hap(ipplot,itor));
                        hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                        name = 'Mesh';
                        title(name), axis equal
                    end
                elseif itor ==1 && ~gatherplot
                    figure(iplot)
                    hold on, plotMesh(Mesh.X,Mesh.T,elemType); axis on, box on
                    name = 'Mesh'; title(name)
                    if printout
                        readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                    end
                end
                if plotCons
                    for ii=1:size(uc,2)
                        uplot = uc(:,ii);
                        if cons_dimensional_plots
                            uplot = uc(:,ii)*simulation_parameters.adimensionalization.reference_values_conservative_variables(ii);
                        end
                        uplot = extractSolutionInAtGivenTheta(uplot,Mesh.T,refEl,refElTor,tpos(itor));
                        if gatherplot
                            icplot = icplot +1; figure(fc(itor)), axes(hac(icplot,itor)); 
                        else
                            iplot = iplot +1; figure(iplot)
                        end
                        hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref,cont,0,0,ccaxis(ii,:));axis off
                        name = simulation_parameters.physics.conservative_variable_names{ii};
                        [num,dem]=rat(tpos(itor)/pi,1e-3); axis equal
                        if num <1e-5
                            title(strcat(name,' - Plane \phi = 0'))
                        else
                            title(strcat(name,' - Plane \phi = ',num2str(num),'\pi','/',num2str(dem)))
                        end
                        if printout
                            readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                        end
                        drawnow
                    end
                end
                if plotPhys
                    for ii=1:size(up,2)
                        uplot = up(:,ii);
                        if phys_dimensional_plots
                            uplot = up(:,ii)*simulation_parameters.adimensionalization.reference_values_physical_variables(ii);
                        end
                        uplot = extractSolutionInAtGivenTheta(uplot,Mesh.T,refEl,refElTor,tpos(itor));
                        if gatherplot
                            ipplot = ipplot +1; figure(fp(itor)), axes(hap(ipplot,itor));
                        else
                            iplot = iplot +1; figure(iplot)
                        end
                        hold on, plotSolution(Mesh.X,Mesh.T,uplot,refEl,nref,cont,0,0,cpaxis(ii,:));axis off
                        name = simulation_parameters.physics.physical_variable_names{ii};
                        [num,dem]=rat(tpos(itor)/pi,1e-3); axis equal
                        if num <1e-5
                            title(strcat(name,' - Plane \phi = 0'))
                        else
                            title(strcat(name,' - Plane \phi = ',num2str(num),'\pi','/',num2str(dem)))
                        end
                        if printout
                            readyforprintjpeg([8 6],16,[],[],1,[],[],path2save,[name,num2str(testnumber)])
                        end
                        drawnow
                    end
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


if (rmp || ripple) && strcmpi(solname(4:5),'3D') %TODO test and subplot
    newiplot = iplot;
    for ip = 1:nproc
        iplot = startPlot + newiplot;
        % Load results
        if parallel
            meshname = [solname(7:pos) '_' num2str(iproc) '_' num2str(nproc) '.h5'];
            solname_comp = [solname '_ip' num2str(iproc) '_it1_np' num2str(nproc) '.h5'];
        else
            meshname = [solname(7:pos), '.h5'];
            solname_comp = [solname, '.h5'];
        end
        HDF5load([solpath,solname_comp]);
        Mesh = HDF5load([meshpath,meshname]);
        Mesh.lscale = simulation_parameters.adimensionalization.length_scale;
        % Number of equations
        Neq = simulation_parameters.Neq;
        %Toroidal angle considered
        theta = simulation_parameters.numerics.Max_extention_in_the_toroidal_direction;
        %tpos = linspace(0+dth,theta-dth,ntpos);
        
        if contains(solname,'west','IgnoreCase',true)
            R0 = 0;
        else
            R0 = simulation_parameters.geometry.Major_radius;
        end
        
        % Correcting mess up with refElem
        if elemType==1
            elemType=0; %quads
        elseif elemType==0
            elemType=1; %triangle
        end
        
        if elemType == 0
            faceNodes = faceNodes_aux_quads(size(Mesh.T,2));
        elseif elemType == 1
            faceNodes = faceNodes_aux(size(Mesh.T,2));
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
        nt = size(magnetic_field,1)/Mesh.Nnodes; % toroidal number of mesh to plot
        smat = size(Mesh.X,1)*nt;
        X3D = zeros(smat,3);
        T3D = zeros(size(Mesh.T,1)*nt,size(Mesh.T,2));
        
        r = sqrt((Mesh.X(:,1) - R0).^2 + Mesh.X(:,2).^2);
        Z = Mesh.X(:,2);
        R = r.*cos(asin(Z./r));
        
        for i=1:nt
            tor_angle = theta*i/nt;
            X3D((i-1)*size(Mesh.X,1)+1:i*size(Mesh.X,1),1) = R.*cos(tor_angle);
            X3D((i-1)*size(Mesh.X,1)+1:i*size(Mesh.X,1),2) = R.*sin(tor_angle);
            X3D((i-1)*size(Mesh.X,1)+1:i*size(Mesh.X,1),3) = Z;
            T3D((i-1)*size(Mesh.T,1)+1:i*size(Mesh.T,1),:) = Mesh.T + (i-1)*size(Mesh.X,1);
        end
        %Connectivity for the faces
        patchFaces = T3D(:,oFaceNodes);
        %Plot 3D mesh
        iplot = iplot +1;
        figure(iplot)
        patchHandle = patch('Faces',patchFaces,'Vertices',X3D,'FaceColor','none','EdgeAlpha',1,'linewidth',0.4, 'linestyle',':');
        
        % Plot of the coils
        % Load coils
        if ripple
            n = numel(coils_ripple(:,1)); % 2 positions per line (3 coordinates)
            coils_ripple = coils_ripple.*Mesh.lscale;
            for i=1:n % segments
                line([coils_ripple(i,1) coils_ripple(i,2)],[coils_ripple(i,3) coils_ripple(i,4)],[coils_ripple(i,5) coils_ripple(i,6)], 'color', 'blue', 'linewidth',2);
            end
        end
        if rmp
            tmp = size(coils_rmp);
            nbRow = tmp(1);
            coils_rmp = coils_rmp.*Mesh.lscale;
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
        hold on, box on
    end
    % stop
    % Plot magnetic fields (perturbation)
    dimB = size(magnetic_field,1);
    ntposB = dimB/Mesh.Nnodes;
    tposB = linspace(0+dth,theta-dth,ntposB);
    
    newiplot = iplot;
    for ip = 1:nproc
        iplot = startPlot + newiplot;
        % Load results
        if parallel
            meshname = [solname(7:pos) '_' num2str(iproc) '_' num2str(nproc) '.h5'];
            solname_comp = [solname '_ip' num2str(iproc) '_it1_np' num2str(nproc) '.h5'];
        else
            meshname = [solname(7:pos), '.h5'];
            solname_comp = [solname, '.h5'];
        end
        HDF5load([solpath,solname_comp]);
        Mesh = HDF5load([meshpath,meshname]);
        Mesh.lscale = simulation_parameters.adimensionalization.length_scale;
        % Number of equations
        Neq = simulation_parameters.Neq;
        %Toroidal angle considered
        theta = simulation_parameters.numerics.Max_extention_in_the_toroidal_direction;
        %tpos = linspace(0+dth,theta-dth,ntpos);
        
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
        refEl = createReferenceElement(elemType,size(Mesh.T,2));
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
        ntposB = dimB/Mesh.Nnodes;
        tposB = linspace(0+dth,theta-dth,ntposB);
        clim1 = 0;
        clim2 = 1e-3;
        
        if average
            BRmean = zeros(Mesh.Nnodes,ntposB);
            BZmean = zeros(Mesh.Nnodes,ntposB);
            BPhimean = zeros(Mesh.Nnodes,ntposB);
            for itor = 1:ntposB
                BRmean(:,itor) = BR((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes);
                BZmean(:,itor) = BZ((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes);
                BPhimean(:,itor) = BPhi((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes);
            end
            BRmean = mean(BRmean,2);
            BZmean = mean(BZmean,2);
            BPhimean = mean(BPhimean,2);
            for itor = 1:ntposB
                %sol = sqrt((BR((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes)-BRmean).^2 + (BZ((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes)-BZmean).^2 + (BPhi((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes)-BPhimean).^2);
                sol = sqrt((BRmean).^2 + (BZmean).^2 + (BPhimean).^2);
                iplot = iplot +1;
                figure(iplot), hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,sol,refEl,nref); axis off, box on
                name = '|B - <B>|';
                [num,dem]=rat(tpos(itor)/pi,1e-3);
                if num < 1e-3
                    title(strcat(name,'- Plane \phi = 0'))
                else
                    title(strcat(name,'- Plane \phi = ',num2str(num),'\pi','/',num2str(dem)))
                end
                caxis([clim1, clim2]);
            end
        else
            for itor = 1:ntposB
                sol = Bmodpert((itor-1)*Mesh.Nnodes+1:itor*Mesh.Nnodes);
                iplot = iplot +1;
                figure(iplot), hold on, plotSolution(Mesh.X/Mesh.lscale,Mesh.T,sol,refEl,nref); axis off, box on
                name = '|Bpert|';
                [num,dem]=rat(tpos(itor)/pi,1e-3);
                if num < 1e-3
                    title(strcat(name,'- Plane \phi = 0'))
                else
                    title(strcat(name,'- Plane \phi = ',num2str(num),'\pi','/',num2str(dem)))
                end
                caxis([clim1, clim2]);
            end
        end
    end
end % rmp/ripple
