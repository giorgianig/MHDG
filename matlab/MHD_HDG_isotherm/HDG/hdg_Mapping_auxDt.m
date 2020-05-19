function [L L0 Lro Lf U U0 Uro Uf P P0 Pro Pf] = hdg_Mapping(...
    flipFace,Nv,nv,Np,A,M,B,BtmLinvA,C,Cv,H,f_conv,D,E,Mp,O,R,W,dt,u0,...
    C_dir,M_dir_prec,O_dir,force)
% matrix multiplication to get the mapping

% number of elements
Ne = size(flipFace,1);

% initialization
L = zeros(4*Nv,6*nv,Ne);
L0 = zeros(4*Nv,1,Ne);
Lro = zeros(4*Nv,1,Ne);
Lf = Lro;
U = zeros(2*Nv,6*nv,Ne);
U0 = zeros(2*Nv,1,Ne);
Uro = zeros(2*Nv,1,Ne);
Uf = Uro;
P = zeros(Np,6*nv,Ne);
P0 = zeros(Np,1,Ne);
Pro = zeros(Np,1,Ne);
Pf = Pro;

% loop in elements
for iElem = 1:Ne
    
    ind_u = (iElem-1)*2*Nv + (1:2*Nv);
    u0e = u0(ind_u);
    
    % elemental matrices
    [Le L0e Lroe Lfe Ue U0e Uroe Ufe Pe P0e Proe Pfe] = elementalMatrices...
        (A(:,:,iElem),M(:,:,iElem),B(:,:,iElem),BtmLinvA(:,:,iElem),...
        C(:,:,iElem),Cv(:,:,iElem),H(:,:,iElem),...
        f_conv(:,:,iElem),D(:,:,iElem),E(:,:,iElem),Mp(:,:,iElem),...
        O(:,:,iElem),R(:,:,iElem),W(:,:,iElem),dt,u0e,...
        C_dir(:,iElem),M_dir_prec(:,iElem),O_dir(:,iElem),force(:,iElem));
    
    % local assembly indexes
    flipFace_e = flipFace(iElem,:);
    ind_1_v_L = (1:2*nv);
    ind_2_v_L = 2*nv + (1:2*nv);
    ind_3_v_L = 4*nv + (1:2*nv);
    
    if flipFace_e(1)
        Le(:,ind_1_v_L) = fliplr2(Le(:,ind_1_v_L));
        Pe(:,ind_1_v_L) = fliplr2(Pe(:,ind_1_v_L));
        Ue(:,ind_1_v_L) = fliplr2(Ue(:,ind_1_v_L));
    end
    if flipFace_e(2)
        Le(:,ind_2_v_L) = fliplr2(Le(:,ind_2_v_L));
        Pe(:,ind_2_v_L) = fliplr2(Pe(:,ind_2_v_L));
        Ue(:,ind_2_v_L) = fliplr2(Ue(:,ind_2_v_L));
    end
    if  flipFace_e(3)
        Le(:,ind_3_v_L) = fliplr2(Le(:,ind_3_v_L));
        Pe(:,ind_3_v_L) = fliplr2(Pe(:,ind_3_v_L));
        Ue(:,ind_3_v_L) = fliplr2(Ue(:,ind_3_v_L));
    end
    
    % store mapping
    L(:,:,iElem) = Le;
    L0(:,:,iElem) = L0e;
    Lro(:,:,iElem) = Lroe;
    Lf(:,:,iElem) = Lfe;
    U(:,:,iElem) = Ue;
    U0(:,:,iElem) = U0e;
    Uro(:,:,iElem) = Uroe;
    Uf(:,:,iElem) = Ufe;
    P(:,:,iElem) = Pe;
    P0(:,:,iElem) = P0e;
    Pro(:,:,iElem) = Proe;
    Pf(:,:,iElem) = Pfe;
end

%% Elemental matrices
function [ LL LL0 LLro LLf UU UU0 UUro UUf PP PP0 PPro PPf] = elementalMatrices...
    (A,M,B,BtmLinvA,C,Cv,H,fconv,D,E,Mp,O,R,W,dt,u0e,C_dir,M_dir_prec,O_dir,force)

% first set
Mu = M + dt*(- Cv + D - BtmLinvA*B); 
Mu_tilde = dt*(E - H - BtmLinvA*C);
Mu0 = M*u0e + dt*fconv;
M_dir = dt*M_dir_prec;
Mp = Mp*dt;

% second set
% invMu = inv(Mu); 
RinvMu = R/Mu;
Gp = RinvMu*Mp;
Gu_tilde = RinvMu*Mu_tilde-O;
Gu = RinvMu*Mu0;
Gf = RinvMu*(M_dir+force*dt)-O_dir;

% third set
K = [Gp W; W' 0];
invK = K\eye(size(K));

% mapping for the pressure
PP = invK(1:end-1,1:end-1)*Gu_tilde;
PP0 = invK(1:end-1,1:end-1)*Gu;
PPro = invK(1:end-1,end);
PPf = invK(1:end-1,1:end-1)*Gf;

% mapping for the velocity
UU = Mu\(Mu_tilde-Mp*PP);
UU0 = Mu\(Mu0-Mp*PP0);
UUro = -Mu\(Mp*PPro);
UUf = Mu\(M_dir+dt*force-Mp*PPf);

% mapping for the velocity gradient
LL = A\(C-B*UU);
LL0 = -A\(B*UU0);
LLro = -A\(B*UUro);
LLf = A\(C_dir-B*UUf);

%% additional routines

function res = fliplr2(A)
% [ 1 1 1      [ 5 5 5 
%   2 2 2        6 6 6
%   3 3 3  ===>  3 3 3
%   4 4 4        4 4 4
%   5 5 5        1 1 1
%   6 6 6 ]      2 2 2 ]
res = zeros(size(A));
res(:,end-1:-2:1) = A(:,1:2:end-1);
res(:,end:-2:2) = A(:,2:2:end);

