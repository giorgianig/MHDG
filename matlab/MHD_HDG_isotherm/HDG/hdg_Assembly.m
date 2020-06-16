function [KK,f] = hdg_Assembly...
    (Fcon,nv,LL,LL0,UU,UU0,Df,Ef,Lf,Hf,Qf,fH,F_dir,nIntFaces)

% mesh data
global neq
nf  = size(Fcon,2);    % number of faces in one element
Ne = size(Fcon,1);                     % number of elements
nOfFaces = max(max(Fcon));          % number of faces
nDirFaces = sum(sum(F_dir));

% allocation and initialization
v_unk = (nOfFaces-nDirFaces)*nv*neq;
dim = v_unk;
allocation = nf^2*Ne*(neq*nv)^2 - nv*nDirFaces + 3*2*2*nv*Ne;

aux_ones_v = ones(1,neq*nf*nv);
index = 0;
I = zeros(allocation,1);
J = zeros(allocation,1);
K = zeros(allocation,1);
f = zeros(dim,1);

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
ind_v_G = ind_v_L;
for iface = 1:nf
    ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

% loop in elements
for iElem = 1:Ne

    Fcone = Fcon(iElem,:);
    aux_ex = Fcone > nIntFaces;
    aux_dir = F_dir(iElem,:);
    aux_ass = true(neq*nf*nv,1);
    for iface = 1:nf
        ind_v_G(iface,:) = (Fcone(iface)-1)*neq*nv + (1:neq*nv);
        aux_ass(ind_v_L(iface,:)) = aux_dir(iface);
    end
%     if aux_ex(1) 
%         Hf(ind_1_v_L,ind_1_v_L,iElem) = 0;
%     end
%     if aux_ex(2)
%         Hf(ind_2_v_L,ind_2_v_L,iElem) = 0;
%     end
%     if aux_ex(3)
%         Hf(ind_3_v_L,ind_3_v_L,iElem) = 0;
%     end
    % global assembly indexes for the velocity
%     ind_1_v_G = (Fcone(1)-1)*2*nv + (1:2*nv);
%     ind_2_v_G = (Fcone(2)-1)*2*nv + (1:2*nv);
%     ind_3_v_G = (Fcone(3)-1)*2*nv + (1:2*nv);
    

    % faces matrices
    KKe_v = Df(:,:,iElem)*UU(:,:,iElem) + Hf(:,:,iElem)-Ef(:,:,iElem) -...
            (Lf(:,:,iElem)-Qf(:,:,iElem))*LL(:,:,iElem);
%     ffe   = Df(:,:,iElem)*U0(:,iElem)+Hfdir(:,iElem);
     ffe   = -Df(:,:,iElem)*UU0(:,iElem) + fH(:,iElem) + ...
            (Lf(:,:,iElem)-Qf(:,:,iElem))*LL0(:,iElem);
% Ke_fort = readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG/test/K_' num2str(iElem) '.txt']);
% fe_fort = readMatTxt(['/home/giorgio/Dropbox/Fortran/MHDG/test/f_' num2str(iElem) '.txt']);
% disp(['Elem' num2str(iElem) ' err K: ' num2str(max(max(abs(Ke_fort-KKe_v))))])
% disp(['Elem' num2str(iElem) ' err f: ' num2str(max(abs(fe_fort-ffe)))])

    % assembly 
    indv_m = col(ind_v_G');
    indv_m = indv_m(~aux_ass);
    indv_transp_m = transpose(indv_m);
    aux_row = indv_m(:,aux_ones_v(~aux_ass));
    aux_col = indv_transp_m(aux_ones_v(~aux_ass),:);
    index = index(end) + (1:numel(indv_m)^2);
    I(index) = aux_row(:);
    J(index) = aux_col(:);
    K(index) = KKe_v(~aux_ass,~aux_ass);
    f(indv_m) = f(indv_m) + ffe(~aux_ass);

end

% check allocation
if size(I,1)>allocation
    error('size overpassed')
end

% create sparse matrix
KK = sparse(I(I~=0),J(I~=0),K(I~=0),dim,dim);
