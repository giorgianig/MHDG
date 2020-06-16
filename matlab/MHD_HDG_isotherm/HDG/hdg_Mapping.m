function [LL,LL0,UU,UU0] = hdg_Mapping(flipFace,Nv,nv,M,G,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
                              B,C,C_dir,invL,P,Pb,Q,Qb)
% matrix multiplication to get the mapping

global neq Mesh

% number of faces per element
nf  = size(flipFace,2);  

% number of dimensions
nd = size(Mesh.X,2);

% number of elements
Ne = size(flipFace,1);

% initialization
LL = zeros(neq*nd*Nv,nf*nv*neq,Ne);
LL0 = zeros(neq*nd*Nv,Ne);
UU = zeros(neq*Nv,neq*nf*nv,Ne);
UU0 = zeros(neq*Nv,Ne);

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
for iface = 1:nf
    ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

% loop in elements
for iElem = 1:Ne
    
    ind_u = (iElem-1)*neq*Nv + (1:neq*Nv);
    u0e = u0(ind_u);
    
    % elemental matrices
    [LLe,LL0e,UUe,UU0e] = elementalMatrices(M(:,:,iElem),G(:,:,iElem),Cv(:,:,iElem),D(:,:,iElem),...
        E(:,:,iElem),Edir(:,iElem),H(:,:,iElem),Hdir(:,iElem),dt,u0e,force(:,iElem),...
        B(:,:,iElem),C(:,:,iElem),C_dir(:,iElem),...
        invL(:,:,iElem),P(:,:,iElem),Pb(:,:,iElem),Q(:,:,iElem),Qb(:,:,iElem));
    
    flipFace_e = flipFace(iElem,:);
    for iface = 1:nf
        if flipFace_e(iface)
            UUe(:,ind_v_L(iface,:)) = UUe(:,ind_v_L(iface,perm));
            LLe(:,ind_v_L(iface,:)) = LLe(:,ind_v_L(iface,perm));
        end
    end
        
    % store mapping
    LL(:,:,iElem) = LLe;
    LL0(:,iElem) = LL0e;    
    UU(:,:,iElem) = UUe;
    UU0(:,iElem) = UU0e;
end

%% Elemental matrices
function [LL,LL0,UU,UU0] = elementalMatrices(M,G,Cv,D,E,Edir,H,Hdir,dt,u0e,f,...
                                     B,C,Cdir,invL,P,Pb,Q,Qb)
PQinvL = ((P-Pb)-(Q-Qb))*invL;

Mu = M/dt-Cv+D+G-PQinvL*B;
Mu_tilde = H-E+PQinvL*C;
Mu0  = f+Edir-Hdir-PQinvL*Cdir+M/dt*u0e;

UU = -Mu\Mu_tilde;
UU0 = Mu\Mu0;
LL = invL*( -B*UU+C);
LL0 = invL*( -B*UU0+Cdir);
