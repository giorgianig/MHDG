% Generate square meshes of high-order quadrangles
close all
clear

%% Define geometry
r0 = 0.5;    % inner radius
re = 1;       % outer radius
rs  = 0.75;  % separatrix radius

%% Define mesh
eltype = 0;    % element type -  0: quads, 1: tria
nonst  = 0;   % alignement with the magnetic field - 0: aligned; 1: non aligned (apply a random displacement to nodes)
delta = 0.0;   % in case on non-alignement: magnitude of the  displacement
ne1 = 3;      % number of layers in the closed-line zone
ne2 = 3;    % number of layers in the open-line zone (SOL)
elong = 3;   % elongation of the elements in the poloidal direction

% Refinement around the limiter
nlref = 1;     % refinement on the limiter: number of original layers in the refinement
lref = 2;       % refinement on the limiter: the number of final layers corresponds to the original one mulptiplyed by this
factorl = 3;   % refinement on the limiter: intensity

% Refinement at the separatrix
nsref = 1;    % refinement on the separatrix: number of original layers in the refinement
sref = 2;      % refinement on the separatrix: the number of final layers corresponds to the original one mulptiplyed by this
factors = 1; % refinement on the separatrix: intensity

%% Define nesting and polynomial degrees
npol = 4;%1:8;
nsiz = 1;%1:7;

%% Working
elemInfo.type = eltype;
if eltype==0
    Nvert = 4;
    elemInfo.nOfNodes = 4;
    elemInfo.faceNodes = [1 2; 2 3; 3 4; 4 1];
else
    Nvert = 3;
    elemInfo.nOfNodes = 3;
    elemInfo.faceNodes = [1 2; 2 3; 3 1];
end
NNodesLin = elemInfo.nOfNodes;
refEl = createReferenceElement(eltype,elemInfo.nOfNodes);
elemInfo.faceNodes1d =  [1 2];
path2save = '/home/giorgio/Matlab/Meshes/Meshes_2D/ConvergenceMeshes/';

if eltype==0
    if nonst
        nameRoot = 'CircLimAlign_Quads_Nonst_';
    else
        nameRoot = 'CircLimAlign_Quads_';
    end
else
    if nonst
        nameRoot = 'CircLimAlign_TriaNONAlign_';
    else
        nameRoot = 'CircLimAlign_TriaAlign_';
    end
end


s =  (rs-r0)/(re-r0);

% number of divisions and element size
% in the poloidal direction withoud refinement
% on the limiter
nx = floor(24*ne2/elong);
hx = 1/nx;

% refined part around the limiter
xs = hx*nlref;
nxref = nlref*lref;

% refined part around the separatrix
hy1 = s/ne1;
ys1 = hy1*nsref;
hy2 = s/ne2;
ys2 = hy2*nsref;
nyref = nsref*sref;

for ic=nsiz
    %     nx = 6*2^(ic-1);
    
    for ip =npol
        %% Create linear mesh
        
        %***** Closed field lines zone************
        % Non refined part
        [X,T]     = CreaMalla(eltype,NNodesLin,xs,1-xs,0    ,s-ys1,nx-2*nlref,ne1-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,xs,1-xs,s-ys1,s    ,nx-2*nlref,nyref);
        q = (s-X3(:,2))/ys1;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s-r*ys1;
        [X,T] = sawMeshes(X,T,X3,T3,s-ys1,2);
        clear X3 T3
        
        % Part refined on the limiter sx
        [X2,T2] = CreaMalla(eltype,NNodesLin,0,xs,0     ,s-ys1,nxref,ne1-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,0,xs,s-ys1,s     ,nxref,nyref);
        q = (s-X3(:,2))/ys1;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s-r*ys1;
        [X2,T2] = sawMeshes(X2,T2,X3,T3,s-ys1,2);
        clear X3 T3
        q = X2(:,1)/xs;
        r = (exp(factorl*q)-1)/(exp(factorl)-1);
        X2(:,1) = r*xs;
        [X,T] = sawMeshes(X,T,X2,T2,xs,1);
        clear X2 T2
        
        % Part refined on the limiter dx
        [X2,T2] = CreaMalla(eltype,NNodesLin,1-xs,1,0      ,s-ys1,nxref,ne1-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,1-xs,1,s-ys1,s       ,nxref,nyref);
        q = (s-X3(:,2))/ys1;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s-r*ys1;
        [X2,T2] = sawMeshes(X2,T2,X3,T3,s-ys1,2);
        clear X3 T3
        q = (1-X2(:,1))/xs;
        r = (exp(factorl*q)-1)./(exp(factorl)-1);
        X2(:,1) = 1-r*xs;
       
        
        [X,T] = sawMeshes(X,T,X2,T2,1-xs,1);
        clear X2 T2
        

        % ******Open field lines zone (SOL)*******
        [X1,T1] = CreaMalla(eltype,NNodesLin,xs,1-xs,s+ys2,1      ,nx-2*nlref,ne2-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,xs,1-xs,s       ,s+ys2,nx-2*nlref,nyref);
        q = (X3(:,2)-s)/ys2;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s+r*ys2;
        [X1,T1] = sawMeshes(X1,T1,X3,T3,s+ys2,2);
        clear X3 T3
        
        % Part refined on the limiter sx
        [X2,T2] = CreaMalla(eltype,NNodesLin,0,xs,s+ys2,1      ,nxref,ne2-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,0,xs,s       ,s+ys2,nxref,nyref);
        q = (X3(:,2)-s)/ys2;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s+r*ys2;
        [X2,T2] = sawMeshes(X2,T2,X3,T3,s+ys2,2);
        clear X3 T3
        q = X2(:,1)/xs;
        r = (exp(factorl*q)-1)/(exp(factorl)-1);
        X2(:,1) = r*xs;
        [X1,T1] = sawMeshes(X1,T1,X2,T2,xs,1);
        clear X2 T2
        
        % Part refined on the limiter dx
        [X2,T2] = CreaMalla(eltype,NNodesLin,1-xs,1,s+ys2,1,nxref,ne2-nsref);
        [X3,T3] = CreaMalla(eltype,NNodesLin,1-xs,1,s       ,s+ys2,nxref,nyref);
        q = (X3(:,2)-s)/ys2;
        r = (exp(factors*q)-1)/(exp(factors)-1);
        X3(:,2) = s+r*ys2;
        [X2,T2] = sawMeshes(X2,T2,X3,T3,s+ys2,2);
        clear X3 T3
        q = (1-X2(:,1))/xs;
        r = (exp(factorl*q)-1)./(exp(factorl)-1);
        X2(:,1) = 1-r*xs;
        [X1,T1] = sawMeshes(X1,T1,X2,T2,1-xs,1);
        clear X2 T2
        
        % Saw the two zones
        [X,T] = sawMeshes(X,T,X1,T1,s,2);
        clear X1 T1
        
        np = size(X,1);
        
        % Generate boundary
        [Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT, elementFaceInfo] =...
            setBoundarySquareMesh(X,T,refEl.faceNodes);
        
        %         figure, plotMesh(X,T,0)
        %         hold on, plot(X(T(47,:),1),X(T(47,:),2),'ro'),
        %         text(X(T(47,:),1),X(T(47,:),2),num2str(transpose(1:size(T,2))),'fontsize',16)
        %         stop
        
        if nonst
            %% Perturbation of the coordinates
            ind_int = setdiff( (1:np)',unique([Tb_UP;Tb_DOWN;Tb_LEFT;Tb_RIGHT]));
            nint = nnz(ind_int);
            %             pert = rand(nint,1);
            pert = col((-1).^(1:nint));
            X(ind_int,1) = X(ind_int,1) + 2*(0.5-pert)*hx*delta;
            X(ind_int,2) = X(ind_int,2) + 2*(0.5-pert)*hy*delta;
        end
        
        
        
        %         figure, plotMesh(X,T)
        
        
        
        if ip>1
            %% Increase order
            [X, T, Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT, ...
                elemInfo, elementFaceInfo] = setOrderLinearMesh...
                (X, T, Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT,...
                elemInfo, elementFaceInfo, ip,eltype);
        end
        
        
        %% Transformation to circular geometry
        X_pol = zeros(size(X));
        X_pol(:,1) = ((re-r0)*X(:,2)+r0).*cos(2*pi*X(:,1));
        X_pol(:,2) = ((re-r0)*X(:,2)+r0).*(-sin(2*pi*X(:,1)));
        
        % Apply 90 rotation in clockwise direction
        X_pol = transpose( [0 1;-1 0]*X_pol');
        
        %% Remove doubled nodes
        unique_Tb_LEFT = unique(Tb_LEFT);
        unique_Tb_RIGHT = unique(Tb_RIGHT);
        [~, sort_LEFT] = sort(X(unique_Tb_LEFT,2));
        [~, sort_RIGHT] = sort(X(unique_Tb_RIGHT,2));
        unique_Tb_LEFT = unique_Tb_LEFT(sort_LEFT);
        unique_Tb_RIGHT = unique_Tb_RIGHT(sort_RIGHT);
        aux_elim = unique_Tb_LEFT > unique_Tb_RIGHT;
        nodes_to_eliminate = zeros(size(unique_Tb_LEFT));
        nodes_to_substitute = nodes_to_eliminate;
        nodes_to_eliminate(aux_elim) = unique_Tb_LEFT(aux_elim) ;
        nodes_to_eliminate(~aux_elim) = unique_Tb_RIGHT(~aux_elim);
        nodes_to_substitute(aux_elim) = unique_Tb_RIGHT(aux_elim);
        nodes_to_substitute(~aux_elim) = unique_Tb_LEFT(~aux_elim);
        T_provv = T;
        Tb_UP_provv = Tb_UP;
        Tb_DOWN_provv = Tb_DOWN;
        for inode = 1:length(nodes_to_substitute)
            node = nodes_to_eliminate(inode);
            T_provv(T==node) = nodes_to_substitute(inode);
            Tb_UP_provv(Tb_UP==node)  = nodes_to_substitute(inode);
            Tb_DOWN_provv(Tb_DOWN==node)  = nodes_to_substitute(inode);
        end
        X_pol(nodes_to_eliminate,:) = [];
        
        % reset T
        nodes = unique(T_provv);
        aux_nodes(nodes) = 1:numel(nodes);
        T = aux_nodes(T_provv);
        Tb_UP = aux_nodes(Tb_UP_provv);
        Tb_DOWN = aux_nodes(Tb_DOWN_provv);
        
        % remove fields
        elementFaceInfo.OUT = elementFaceInfo.UP;
        elementFaceInfo.IN = elementFaceInfo.DOWN;
        elementFaceInfo = rmfield(elementFaceInfo, {'LEFT' 'RIGHT','UP','DOWN'});
        Tb_OUT = Tb_UP;
        Tb_IN = Tb_DOWN;
        clear X Tb_LEFT Tb_RIGHT Tb_UP_provv Tb_DOWN_provv Tb_DOWN Tb_UP...
            T_provv unique_Tb_LEFT unique_Tb_RIGHT
        X = X_pol;
        
        
        %% Add limiter
        tol = 1e-10;
        % Find nodes that lay on the limiter
        limnodes = find(all([X(:,1)<tol, X(:,1)>-tol, X(:,2)<=-rs],2));
        limvert = limnodes(all([X(limnodes,2)<-rs+tol,X(limnodes,2)>-rs-tol],2));
        limnodes_rep = setdiff(limnodes,limvert); % nodes to replicate
        
        % Find elements connected to the limiter
        limels = find(sum(ismember(T(:,1:Nvert),limnodes),2)==2);
        % Elements that lay to the right of the limiter
        limels_r = limels(all(transpose(reshape(X(T(limels,1:Nvert)',1),[Nvert,numel(limels)]))>=-tol,2));
        % Elements that lay on the left of the limiter
        limels_l = setdiff(limels,limels_r);
        
        % Add nodes
        X = [X; X(limnodes_rep,:)];
        limnodes_in = max(max(T))+(1:numel(limnodes_rep));
        
        
        Tb_OUT_old = Tb_OUT;        
        %
        elementFaceInfo.LIM = zeros(numel(limels),2);
        Tb_LIM = zeros(numel(limels),numel(elemInfo.faceNodes1d));
        for i = 1:numel(limels_l)
            iel = limels_l(i);
            nd = ismember(T(iel,1:Nvert),limnodes);
            face = findface(nd);
            Tb_LIM(i,:) = T(iel,elemInfo.faceNodes(face,:));
            elementFaceInfo.LIM(i,:) = [iel,face];
        end
        for i = 1:numel(limels_r)
            iel = limels_r(i);
            nd = ismember(T(iel,1:Nvert),limnodes);
            face = findface(nd);
            ndi = find(ismember(T(iel,:),limnodes));
            for ino = 1:numel(ndi)
                if T(iel,ndi(ino))~=limvert
                    if ismember(T(iel,ndi(ino)),Tb_OUT(:,[1,end]))
                        f_out = find(all(ismember(Tb_OUT(:,[1,end]),T(iel,1:Nvert)),2));
                        Tb_OUT(f_out,Tb_OUT(f_out,[1,end])==T(iel,ndi(ino)))=limnodes_in(ismember(limnodes_rep,T(iel,ndi(ino))));
                    end
                    T(iel,ndi(ino)) = limnodes_in(ismember(limnodes_rep,T(iel,ndi(ino))));
                end
            end
            Tb_LIM(i+numel(limels_l),:) = T(iel,elemInfo.faceNodes(face,:));
            elementFaceInfo.LIM(i+numel(limels_l),:) = [iel,face];
        end
        
        
        
        
        
        %% Plot mesh
        if eltype==1
            figure, plotMesh(X,T)
        else
            figure, plotMesh(X,T,0)
        end
        hold on
        plot(X(Tb_LIM,1),X(Tb_LIM,2),'ro')
        
        disp(['Generated a mesh with ', num2str(size(T,1)), ' elements'])
        %% Save mesh
%         save([path2save nameRoot 'Nel' num2str(size(T,1)) '_P' num2str(ip) '.mat'],'X','T','Tb_IN','Tb_OUT','Tb_LIM',...
%             'elementFaceInfo','elemInfo')
        
    end
end