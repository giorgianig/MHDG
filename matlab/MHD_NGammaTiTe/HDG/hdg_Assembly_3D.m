function [KK,f] = hdg_Assembly_3D...
    (refEl,F,F_dir,LL,LL0,UU,UU0,Df,Ef,Lf,Hf,Qf,fH, TUhf,TQhf,Tfhf)

% mesh data
global  ntor neq refElTor
Np2d      = size(refEl.NodesCoord,1);              % Number of points for 2d element
Np1dPol   = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor   = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
N2d       = size(F,1);                                                % Number of 2d elements
nf        = size(refEl.faceNodes,1);                      % Number of faces in the 2d element
Nfl       = Np1dPol*Np1dTor;                                    % Number of nodes for lateral faces
Nfp       = Np2d*2+nf*Np1dPol*Np1dTor;                  % Number of nodes in all the 3d faces of 1 element
Nf        = max(max(F));                                              % Number of faces in the 2d plane
nDirFaces = sum(sum(F_dir));
Ne        = N2d*ntor;                                               % Number of 3d elements

% Number of lateral unknown faces
nunk_face_tor  = (Nf-nDirFaces)*ntor;

% Number of unknown faces at extremities
nunk_face_pol = N2d*ntor;

% Matrix dimensions
dim = neq*Nfl*nunk_face_tor + neq*Np2d*nunk_face_pol;

% Allocation for assembly vectors
allocation = Ne*(Nfp*neq)^2-(nDirFaces*ntor)*(Nfl*neq)^2;

% Vector initialization
aux_ones = ones(1,Nfp*neq);
index = 0;
I = zeros(allocation,1);
J = zeros(allocation,1);
K = zeros(allocation,1);
f = zeros(dim,1);
ind_ass = zeros(Nfp*neq,1);
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
aux_dir = false(1,nf+2);
% loop in elements
for itor = 1:ntor
    
    %% Loop in 2d elements
    for iel = 1:N2d
        
        %% Some indices
        iElem = (itor-1)*N2d+iel; % 3d numbering of the element
        Fe = F(iel,:);
        aux_ass = true(size(ind_ass));
        aux_dir(2:end-1) = F_dir(iel,:);
        
        % faces matrices
        KKe_v = Df(:,:,iElem)*UU(:,:,iElem) + Hf(:,:,iElem)-TUhf(:,:,iElem)-Ef(:,:,iElem) -...
                (Lf(:,:,iElem)+TQhf(:,:,iElem)-Qf(:,:,iElem))*LL(:,:,iElem);
         ffe   = -Df(:,:,iElem)*UU0(:,iElem) + fH(:,iElem) -Tfhf(:,iElem) +...
                (Lf(:,:,iElem)+TQhf(:,:,iElem)-Qf(:,:,iElem))*LL0(:,iElem);
      
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*[(iel-1)*Np2d,N2d*Np2d+...
            (Fe-1)*Nfl,N2d*Np2d+(Nf-nDirFaces)*Nfl+(iel-1)*Np2d];
        if itor==ntor
            delta(end) = 1+(iel-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            ind_ass(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
            aux_ass(ind_loc) = ~aux_dir(iface);
        end
        % assembly
        ind_ass = ind_ass(aux_ass);
        ind_ass_t = transpose(ind_ass);
        aux_row = ind_ass(:,aux_ones(aux_ass));
        aux_col = ind_ass_t(aux_ones(aux_ass),:);
        index = index(end) + (1:numel(ind_ass)^2);
        I(index) = aux_row(:);
        J(index) = aux_col(:);
        K(index) = KKe_v(aux_ass,aux_ass);
        auxffe = ffe(aux_ass);
        for i=1:numel(ind_ass)
        f(ind_ass(i)) = f(ind_ass(i)) + auxffe(i);
        end
% f_fort = readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/fe' num2str(iElem) '.txt']);
% max(abs(f-f_fort))
        
    end
end
% check allocation
if size(I,1)>allocation
    error('size overpassed')
end

% create sparse matrix
KK = sparse(I(I~=0),J(I~=0),K(I~=0),dim,dim);
