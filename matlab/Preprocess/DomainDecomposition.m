clear
clc
close all
mark = {'bo','ko','go','mo','yo','bo','bo','ko','go','mo','yo','bo'};
% Mesh to load
path2mesh ='../Meshes/';
% meshName  = 'Circle_ONION_3_P10.mat';
% meshName  = 'West_Hugo_h0.02_refCorn0.001_refSep0.01_YesHole_P4.mat';
meshName  = 'CircLimAlign_Quads_Nel1536_P3.mat';
% meshName  = 'CircLimAlign_Quads_Nel128_P3.mat';
% path2mesh = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/';
% meshName  = 'mesh1_P1.mat';
% Number of divisions
ndiv = 12;
divType = 'rot';
% divType = 'x';
% divType = 'y';
theta0 = 1*pi/180;
elemType = 0;

% theta0 = pi/180;
% forbidden values for theta
theta_forb = [270];

if elemType==1 % Quadrangles
    Nvert =3;
elseif elemType==0 % Triangles
    Nvert = 4;
else
    error('Wrong element type')
end
if elemType==0
    elemTypeFort = 1;
elseif elemType==1
    elemTypeFort =0;
end

%% Loading original mesh
load([path2mesh meshName]);
NelTot = size(T,1);
refEl = createReferenceElement(elemType,size(T,2),[]);
aux = whos('Tb_*');
boundaryNames = {aux.name};
nboundaries = numel(boundaryNames);
nv = size(refEl.faceNodes1d,2);
if exist('intfaces','var')
    intFaces_Glob=intfaces;
    extFaces_Glob=extfaces;
    clear intfaces extfaces
else
    [intFaces_Glob,extFaces_Glob] = GetFaces(T(:,1:Nvert),refEl.faceNodes);
    intfaces = intFaces_Glob;
    extfaces = extFaces_Glob;
    save([path2mesh meshName],'intfaces','extfaces','-append')
    clear intfaces extfaces
end
Nf = size(intFaces_Glob,1)+size(extFaces_Glob,1);
Nintfaces = size(intFaces_Glob,1);
mnx = 0.5*(max(X(:,1))+min(X(:,1)));
mny = 0.5*(max(X(:,2))+min(X(:,2)));


%% SCALING to mach the right dimensions for GOLDSTONE
% X = X/0.75*0.287;

figure, plotMesh(X,T,elemType),hold on

%% Initialization
TT = cell(ndiv,1);
XX = TT;
Loc2Glob_no  = TT;
Loc2Glob_el  = TT;
Loc2Glob_fa  = TT;
ReptFaces    = TT;
intfaces_loc = TT;
extfaces_loc = TT;
Tb_loc       = TT;
mapextfa_loc = TT;
flipFaces_GHOST = TT;
GhostFaces    = TT;
GhostElems    = TT;
GhostPro      = TT;
GhostLoc      = TT;
GhelsPro      = TT;
GhelsLoc      = TT;

s1 = 0;
for ib = 1:nboundaries
    name = boundaryNames{ib};
    eval(['TTb_' name(4:end) '=TT;'])
    eval(['s1 = s1+size(Tb_' name(4:end) ',1);']);
end
TTb_CON = TT;
numberOfFaces=zeros(ndiv,1);
nintf=zeros(ndiv,1);


%% Global boundary connectivity TTb
% TTb for loc2glob faces
TTb = zeros(s1,nv);
s1 = 1;
for ibound = 1:nboundaries
    name = boundaryNames{ibound};
    eval(['aux_bound = Tb_' name(4:end) ';'])
    TTb(s1:size(aux_bound,1)+s1-1,:) = aux_bound;
    s1=s1+size(aux_bound,1);
end


%% ***************** STARTING DIVISION ******************************
if strcmp(divType,'rot')
    theta = linspace(0,2*pi,ndiv+1);
    
    % elements baricentre
    x = X(:,1)-mnx; y = X(:,2)-mny;
    % rotation
    rot = [cos(theta0) -sin(theta0);sin(theta0) cos(theta0)];
    aux = rot*transpose([x,y]);
    x = aux(1,:)';y = aux(2,:)';
    
elseif strcmp(divType,'x')
    theta = linspace(min(X(:,1)),max(X(:,1)),ndiv+1);
    % elements baricentre
    x = X(:,1); y = X(:,2);
    
elseif strcmp(divType,'y')
    theta = linspace(min(X(:,2)),max(X(:,2)),ndiv+1);
    % elements baricentre
    x = X(:,1); y = X(:,2);
    
end

%
xb = 1/Nvert*sum(x(T(:,1:Nvert)),2);
yb = 1/Nvert*sum(y(T(:,1:Nvert)),2);

if strcmp(divType,'rot')
    [theta_b,~] = cart2pol(xb,yb);
    theta_b(theta_b<0) = theta_b(theta_b<0)+2*pi;
    % Rotate back
    x = X(:,1)-mnx; y = X(:,2)-mny;
    xb = 1/Nvert*sum(x(T(:,1:Nvert)),2);
    yb = 1/Nvert*sum(y(T(:,1:Nvert)),2);
elseif strcmp(divType,'x')
    theta_b = xb;
elseif strcmp(divType,'y')
    theta_b = yb;
end




%% Balancing divisions

nel_mean = NelTot/ndiv;
disp('Balancing divisions')
nel = zeros(ndiv,1);
thetalim0=theta(2:end-1);
thetalim=thetalim0;
redo = 1;
imax = 1000;
while redo
    redo=0;
    for it=1:ndiv-1
        disp(['Division',num2str(it)])
        theta0 = theta(it);
        theta1 = theta(it+2);
        for i=1:imax
            % compute number of elements current division
            ind_it = all([theta_b>theta(it) theta_b<theta(it+1)],2);
            nel = nnz(ind_it);
            if (nel-nel_mean)>0
                theta1=theta(it+1);
                theta(it+1) = 0.5*(theta0+theta1);
            elseif (nel-nel_mean)<0
                theta0 = theta(it+1);
                theta(it+1) = 0.5*(theta0+theta1);
            else
                disp('Ok')
                break
            end
        end
        if i==imax
            disp('Optimal division not found, retrying')
            redo = 1;
        end
    end
    if any(theta-theta_forb)<5e-3
        theta = theta+0.1;
        redo=1;
    end

end
for it=1:ndiv
    ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
    nel(it) = nnz(ind);
end

%
% for iter=1:10
%     [balance0,nel] = compute_balance(thetalim,theta_b);
%     balance0
%     nel
%     if balance0<100,break,end
%     jac=zeros(1,ndiv);
%     hes=jac;
%     for it=1:ndiv-1
%         delta=1;
%         thetavarp=thetalim;
%         thetavarm=thetalim;
%         for i=1:100
%             dtheta = delta*(theta(it+2)-theta(it+1));
%             thetavarp(it)=thetalim(it)+dtheta;
%             thetavarm(it)=thetalim(it)-dtheta;
%             balanceit = compute_balance(thetavarp,theta_b);
%             balancemt = compute_balance(thetavarm,theta_b);
%             if (balanceit~=balance0)
%                 delta=delta*0.8;
%                 balanceit_ok=balanceit;
%                 balancemt_ok=balancemt;
%                 dtheta_ok=dtheta;
%             else
%                 jac(it+1)=(balanceit_ok-balance0)/dtheta_ok;
%                 hes(it+1)=(balanceit_ok-2*balance0+balancemt_ok)/dtheta_ok^2;
%                  break
%             end
%         end
%     end
%     varth=jac./hes;
%     thetalim=thetalim-varth(2:end);
% end
% stop

% f=@(thetalim)balance(thetalim,theta_b);
% theta_balance = fminsearch(f,thetalim);
% theta = [0,theta_balance,2*pi];
% for it=1:ndiv
%     ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
%     nel(it) = nnz(ind);
% end
% bal = (nel-nel_mean)/nel_mean;
% disp(['Done balancing: balancing achieved:',num2str(bal)])




% plot(xb,yb,'*')
%  nel_mean = NelTot/ndiv;
% disp('Balancing divisions')
% nel = zeros(ndiv,1);
% nel_mean = NelTot/ndiv;
% thetalim=theta(2:end-1);
% f=@(thetalim)balance(thetalim,theta_b);
% theta_balance = fminsearch(f,thetalim);
% theta = [0,theta_balance,2*pi];
% for it=1:ndiv
%     ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
%     nel(it) = nnz(ind);
% end
% bal = (nel-nel_mean)/nel_mean;
% disp(['Done balancing: balancing achieved:',num2str(bal)])



% tol = 0.01;
% nel = zeros(ndiv,1);
% nel_mean = NelTot/ndiv;
% for itry =1:100
% %     for it=1:ndiv
% %         ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
% %         nel(it) = nnz(ind);
% %     end
% %     nel
% %     bal = (nel-nel_mean)/nel_mean;
%     [bal,nel] = compute_balance(theta,theta_b);
%     nel
%     if itry==1
%         disp(['Initial balancing: ' num2str(100*bal') ' %'])
%     end
%     if any(abs(bal)>tol)
%         for it = 1:ndiv-1
%             check=true;
%             while check
%                 theta(it+1) = theta(it+1) - bal(it)*0.1*(theta(it+1)-theta(it));
%                 [bal_it,nel_it] = compute_balance_it(theta,theta_b,it,nel_mean);
%                 if (bal(it)~=bal_it)
%                     check=false;
%                 end
%             end
%             [bal,nel] = compute_balance(theta,theta_b);
%         end
%         for ifor = 1:numel(theta_forb)
%             k = abs(theta-theta_forb(ifor))<5e-3;
%             theta(k) = theta(k)+5e-3;
%         end
%     else
%         break
%     end
% %     for it=1:ndiv
% %         ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
% %         nel(it) = nnz(ind);
% %     end
% %     nel
% %     bal = (nel-nel_mean)/nel_mean;
%     if any(abs(bal)>tol)
%         for it = ndiv:-1:2
% %             theta(it) = theta(it) + bal(it)*0.2*(theta(it+1)-theta(it));
%             while check
%                 theta(it+1) = theta(it+1) - bal(it)*0.1*(theta(it+1)-theta(it));
%                 [bal_it,nel_it] = compute_balance_it(theta,theta_b,it,nel_mean);
%                 if (bal(it)~=bal_it)
%                     check=false;
%                 end
%             end
%             [bal,nel] = compute_balance(theta,theta_b);
%         end
%         for ifor = 1:numel(theta_forb)
%             k = abs(theta-theta_forb(ifor))<5e-3;
%             theta(k) = theta(k)+5e-3;
%         end
%     else
%         break
%     end
% end

for it = 1:ndiv
    disp(['Division',num2str(it)])
    %% Make a first division of the mesh without overlapping
    ind = all([theta_b>theta(it) theta_b<theta(it+1)],2);
    loc2glob_el = find(ind);
    
    %% Add ghost elements to the divisions (overlapping)
    disp('Add ghost elements to the division')
    auxgh = [ismember(intFaces_Glob(:,1),loc2glob_el),ismember(intFaces_Glob(:,3),loc2glob_el)];
    indgh = sum(auxgh,2)==1;
    ghel = unique(intFaces_Glob(indgh,[1,3]));
    innext = find(theta_b(ghel)>theta(it+1));
    ind(ghel(innext)) = true;
    loc2glob_el = find(ind);
    GhostElems{it} =double(ismember(loc2glob_el,ghel(innext)));
    
    % plot
    %     if strcmp(divType,'rot')
    %         plot(xb(ind)+mnx,yb(ind)+mny, mark{it})
    %     elseif strcmp(divType,'x')
    %         plot(xb(ind),yb(ind), mark{it})
    %     elseif strcmp(divType,'y')
    %         plot(xb(ind),yb(ind), mark{it})
    %     end
    %
    
    %% Generate X and T for the current division
    disp('Generate X and T for the current division')
    Tprov = T(ind,:);
    nodes = unique(Tprov);
    Xprov = X(nodes,:);
    %     Tfix = Tprov;
    count = 0;
    mapTfix = unique(Tprov);
    invmap = zeros(mapTfix(end),1);
    invmap(mapTfix) = (1:size(Xprov,1))';
    Tfix = invmap(Tprov);
    %     for i = 1:size(Xprov,1)
    %         if mapTfix(i)==i,continue,end
    %         Tfix(Tfix==mapTfix(i)) = i;
    % %         while isempty(find(Tfix(:)==i, 1))
    % %             Tfix(Tfix>i) = Tfix(Tfix>i)-1;
    % %             count = count+1;
    % %         end
    %
    %     end
    nodesfix = unique(Tfix);
    % check that everything is ok so far...
    if ~isequal(nodesfix,transpose(1:max(max(Tfix))))
        error('Problem in fixed connectivity')
    end
    
    
    
    %% Create loc2glob_fa: interior faces
    disp('Create loc2glob_fa: interior faces')
    [intfaces,extfaces] = GetFaces(Tfix(:,1:Nvert),refEl.faceNodes);
    intfaces_loc{it} = intfaces;
    extfaces_loc{it} = extfaces;
    numberOfFaces(it) = size(intfaces,1)+size(extfaces,1);
    nintf(it) = size(intfaces,1);
    loc2glob_fa = zeros(numberOfFaces(it),1);
    reptFaces=loc2glob_fa;
    ghostFaces=loc2glob_fa;
    for j=1:size(intfaces,1)
        faceloc = intfaces(j,:);
        faceglo = [loc2glob_el(faceloc(1)),faceloc(2),...
            loc2glob_el(faceloc(3)),faceloc(4),faceloc(5)];
        %         [a,b] = ismember(faceglo,intFaces_Glob,'rows');
        b = find(sum(abs(intFaces_Glob-repmat(faceglo,[size(intFaces_Glob,1),1])),2)==0);
        %         if ~a,error('Problem in the creation of loc2glob_fa'),end
        if isempty(b),error('Problem in the creation of loc2glob_fa'),end
        loc2glob_fa(j) = b;
    end
    
    
    
    
    %% Generate boundary connectivity for the current division
    disp('Generate boundary connectivity for the current division')
    clg = 1;
    for ib = 1:nboundaries
        name = boundaryNames{ib};
        eval(['Tb = Tb_' name(4:end) ';'])
        Tbprov = zeros(size(Tb));
        count = 1;
        for ifa = 1:size(Tb,1)
            n1 = Tb(ifa,1);
            n2 = Tb(ifa,end);
            if all([ismember(n1,Tprov) ismember(n2,Tprov)])
                Tbprov(count,:) = Tb(ifa,:);
                count = count+1;
                [a,b] = ismember(Tb(ifa,:),TTb,'rows');
                if ~a,error('Problem in the creation of loc2glob_fa'),end
                loc2glob_fa(clg+nintf(it)) = b+Nintfaces;
                clg=clg+1;
            end
        end
        Tbprov = Tbprov(1:count-1,:);
        for nod=col(unique(Tbprov))'
            rep = find(nodes==nod);
            if isempty(rep), error('Something wrong!'), end
            Tbprov(Tbprov==nod) = rep;
        end
        eval(['TTb_' name(4:end) '{it}=Tbprov;'])
    end
    
    
    %% Store subdivision
    Loc2Glob_fa{it} = loc2glob_fa;
    ReptFaces{it} = reptFaces;
    TT{it} = Tfix;
    XX{it} = Xprov;
    Loc2Glob_no{it} = nodes;
    Loc2Glob_el{it} = loc2glob_el;
    GhostFaces{it}    = ghostFaces;
end



% stop
%% Boundary connectivity for the overlapping faces
disp('Boundary connectivity for the overlapping faces')
faceNodes = refEl.faceNodes;
if elemType==0
    faceslin = [1 2; 2 3; 3 4; 4 1];
elseif elemType==1
    faceslin = [1 2; 2 3; 3 1];
end

% if elemType==1 % Quadrangles
%     faceslin = [1 2; 2 3; 3 4; 4 1];
% elseif elemType==0 % Triangles
%     faceslin = [1 2; 2 3; 3 1];
% end
for i=1:ndiv
    disp(['Division',num2str(i)])
    Tp = TT{i};
    loc2glob_fa = Loc2Glob_fa{i};
    loc2glob_el = Loc2Glob_el{i};
    ghostFaces  = GhostFaces{i} ;
    T_CON = zeros(size(Tp,1),size(faceNodes,2));
    flipFace = zeros(size(Tp,1),1);
    extFaces = extfaces_loc{i};
    count = 1;
    clg = find(loc2glob_fa==0,1);
    mapextf = zeros(size(extFaces,1),1);
    for iFa = 1:size(extFaces,1)
        iel = extFaces(iFa,1);
        ifa = extFaces(iFa,2);
        n1 = Tp(iel,faceslin(ifa,1));
        n2 = Tp(iel,faceslin(ifa,2));
        check = 1;
        s = 0;
        for ib = 1:nboundaries
            name = boundaryNames{ib};
            eval(['Tb = TTb_' name(4:end) '{i};'])
            if all([ismember(n1,Tb) ismember(n2,Tb)])
                [a1,b1]=ismember([n1 n2],Tb(:,[1,end]),'rows');
                [a2,b2]=ismember([n2 n1],Tb(:,[1,end]),'rows');
                if a1
                    mapextf(iFa)= b1+s;
                elseif a2
                    mapextf(iFa)= b2+s;
                else
                    error('strange')
                end
                check=0;
                break
            end
            s = s+size(Tb,1);
        end
        if check
            T_CON(count,:) = Tp(iel,faceNodes(ifa,:));
            mapextf(iFa) = count+s;
            ghostFaces(count+s+nintf(i)) = 1;
            [a1,b1]=ismember([loc2glob_el(iel),ifa],intFaces_Glob(:,1:2),'rows');
            [a2,b2]=ismember([loc2glob_el(iel),ifa],intFaces_Glob(:,3:4),'rows');
            if a1
                loc2glob_fa(clg)=b1;
                if a2,error('That is strange'),end
                if loc2glob_el(iel)>intFaces_Glob(b1,3),flipFace(count)=1;end
            elseif a2
                loc2glob_fa(clg)=b2;
                if loc2glob_el(iel)>intFaces_Glob(b2,1),flipFace(count)=1;end
            else
                error('That is very strange')
            end
            count = count+1;
            clg=clg+1;
        end
    end
    % little check
    if any(mapextf==0),error('error in mapextf'),end
    if ~isequal(unique(mapextf)',1:max(mapextf)),error('error in mapextf'),end
    TTb_CON{i}       = T_CON(1:count-1,:);
    flipFaces_GHOST{i} = flipFace(1:count-1);
    Loc2Glob_fa{i}   = loc2glob_fa;
    mapextfa_loc{i}  = mapextf;
    GhostFaces{i}    = ghostFaces;
end

%% Find multiple occurrence of the same face across various proc
disp('Find multiple occurrence of the same face across various proc')
checkFaces = zeros(Nf,1);
for i = 1:ndiv
    disp(['Division',num2str(i)])
    Tp = TT{i};
    loc2glob_el = Loc2Glob_el{i};
    reptFaces   = ReptFaces{i};
    ghostFaces  = GhostFaces{i};
    [intFaces,extFaces] = GetFaces(Tp(:,1:Nvert),refEl.faceNodes);
    
    % Interior faces
    [ckint,pos] = ismember(intFaces_Glob,[loc2glob_el(intFaces(:,1)),intFaces(:,2),...
        loc2glob_el(intFaces(:,3)),intFaces(:,4),...
        intFaces(:,5)],'rows');
    pos = pos(pos~=0);
    reptFaces(pos) = checkFaces(ckint);
    checkFaces(ckint) = 1;
    
    % Exterior faces
    [ckext,pos] = ismember(extFaces_Glob,[loc2glob_el(extFaces(:,1)),extFaces(:,2)],'rows');
    pos = pos(pos~=0);
    reptFaces(size(intFaces,1)+pos) = checkFaces(size(intFaces_Glob,1)+find(ckext));
    checkFaces(size(intFaces_Glob,1)+find(ckext)) = 1;
    
    % Change order in exterior faces
    auxext = reptFaces(size(intFaces,1)+1:end);
    auxext(mapextfa_loc{i} ) = reptFaces(size(intFaces,1)+1:end);
    reptFaces(size(intFaces,1)+1:end) = auxext;
    ghostFaces(reptFaces==1) = 1;
    
    ReptFaces{i} = reptFaces;
    GhostFaces{i}= ghostFaces;
    
    % fix flipFaces dimension
    auxflp = zeros(numel(find(ghostFaces)),1);
    nfl = numel(flipFaces_GHOST{i});
    auxflp(end-nfl+1:end) = flipFaces_GHOST{i};
    flipFaces_GHOST{i} = auxflp;
    
end
clear T_CON Xp Tp Tb* Tfix nodes nodesfix theta theta_b
clear rho_b rep X T Xprov x y xb yb Tb Tprov
clear i ib iel ifa iFa intFaces it loc2glob nod faceNodes
clear faceslin check extFaces n1 n2 aux_bound TTb s1 intfaces_glob


%% Create structure for the communications
disp('Create structure for the communications')
for i = 1:ndiv
    disp(['Division',num2str(i)])
    ghost = find(GhostFaces{i});
    ghost_glob = Loc2Glob_fa{i}(ghost);
    ghostpro  = zeros(size(ghost));
    ghostloc  = ghostpro;
    for j=1:ndiv
        
        if j==i,continue,end
        [a,b] = ismember(ghost_glob,Loc2Glob_fa{j});
        ghostpro(a) = j;
        ghostloc(a)  = b(a);
    end
    GhostPro{i} = ghostpro;
    GhostLoc{i} = ghostloc;
    
    ghels = find(GhostElems{i});
    ghels_glob = Loc2Glob_el{i}(ghels);
    ghelspro  = zeros(size(ghels));
    ghelsloc  = ghelspro;
    for j=1:ndiv
        
        if j==i,continue,end
        [a,b] = ismember(ghels_glob,Loc2Glob_el{j});
        ghelspro(a) = j;
        ghelsloc(a)  = b(a);
    end
    GhelsPro{i} = ghelspro;
    GhelsLoc{i} = ghelsloc;
end

% Draw
% for i=1:ndiv
%     
%     figure,title(['Process ',num2str(i)])
%     plotMesh(XX{i},TT{i},elemType)
%     for ib = 1:nboundaries
%         name = boundaryNames{ib};
%         eval(['Tb = TTb_' name(4:end) '{i};'])
%         hold on
%         plot(XX{i}(Tb,1),XX{i}(Tb,2), mark{ib})
%     end
%     %     Tb = TTb_CON{i};
%     %     hold on
%     %     plot(XX{i}(Tb,1),XX{i}(Tb,2),  'ro')
%     
%     % draw ghost faces
%     ghost = find(GhostFaces{i});
%     nint = size(intfaces_loc{i},1);
%     extfacesperm(mapextfa_loc{i},:) = extfaces_loc{i};
%     for ig=1:numel(ghost)
%         
%         gf = ghost(ig);
%         if gf>nint
%             iel = extfacesperm(gf-nint,1);
%             ifa = extfacesperm(gf-nint,2);
%         else
%             iel = intfaces_loc{i}(gf,1);
%             ifa = intfaces_loc{i}(gf,2);
%         end
%         hold on
%         plot(XX{i}(TT{i}(iel,refEl.faceNodes(ifa,:)),1),XX{i}(TT{i}(iel,refEl.faceNodes(ifa,:)),2),'r*')
%     end
%     
%     
% end

GenerateHDF5_meshes_Parallel
% test