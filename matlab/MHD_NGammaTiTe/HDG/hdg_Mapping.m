function [LL,LL0,UU,UU0] = hdg_Mapping(flipFace,Nv,nv,M,G,Cv,H,Hdir,D,E,Edir,dt,u0,force,...
                              B,C,C_dir,invL,P,Pb,Q,Qb, TU,TUh,TQ,TQh,Tf,Tfh,TUhdir)
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
        invL(:,:,iElem),P(:,:,iElem),Pb(:,:,iElem),Q(:,:,iElem),Qb(:,:,iElem),...
        TU(:,:,iElem),TUh(:,:,iElem),TQ(:,:,iElem),TQh(:,:,iElem),Tf(:,iElem),Tfh(:,iElem),TUhdir(:,iElem));
    
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
                                     B,C,Cdir,invL,P,Pb,Q,Qb,TU,TUh,TQ,TQh,Tf,Tfh,TUhdir)

                                 
% C_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/C.txt');
% E_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/E.txt');
% H_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/H.txt');
% Cv_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/Cv.txt');
% Cv_bef_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/Cv_bef.txt');
% PinvL= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/PinvL.txt');
% MdtCvD= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/MdtCvD.txt');
% MdtCvDG= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/MdtCvDG.txt');
% P_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/P.txt');
% TU_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/TU.txt');
% PTQ_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/PTQ.txt');
% 
% 
% disp(['Error on P:' num2str(max(abs(col((P-Pb)-(Q-Qb))-P_fort(:))))])
% disp(['Error on Cv:' num2str(max(abs(col(Cv)-Cv_bef_fort(:))))])
% disp(['Error on TU:' num2str(max(abs(col(TU)-TU_fort(:))))])


% Energy stuff
P = P+TQ-TQh;
Cv = Cv-TU;
H = H-TUh;
Edir = Edir+TUhdir-Tfh+Tf;

% Old terms without energy
PQinvL = ((P-Pb)-(Q-Qb))*invL;
Mu = M/dt-Cv+D+G-PQinvL*B;
Mu_tilde = H-E+PQinvL*C;
Mu0  = f+Edir-Hdir-PQinvL*Cdir+M/dt*u0e;


% disp(['Error on P+TQ: ' num2str(max(abs(col((P-Pb)-(Q-Qb))-PTQ_fort(:))))])

% disp(['Error on C: ' num2str(max(abs(C(:)-C_fort(:))))])
% disp(['Error on E: ' num2str(max(abs(E(:)-E_fort(:))))])
% disp(['Error on H: ' num2str(max(abs(H(:)-H_fort(:))))])
% disp(['Error on Cv: ' num2str(max(abs(Cv(:)-Cv_fort(:))))])
% disp(['Error on PQinvL:' num2str(max(abs(PQinvL(:)-PinvL(:))))])
% disp(['Error on MdtCvD: ' num2str(max(abs(col(MdtCvD)-col(M/dt-Cv+D))))])
% disp(['Error on MdtCvDG: ' num2str(max(abs(col(MdtCvDG)-col(M/dt-Cv+D+G))))])

% stop
if cond(Mu)>1e14, error('Error: malconditioned elemental matrix'),end
UU = -Mu\Mu_tilde;
UU0 = Mu\Mu0;
LL = invL*( -B*UU+C);
LL0 = invL*( -B*UU0+Cdir);


% Mu_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/Mu.txt');
% Mt_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/Mt.txt');
% % 
% disp(['Error on Mu: ' num2str(max(abs(Mu(:)-Mu_fort(:))))])
% disp(['Error on Mu_tilde: ' num2str(max(abs(Mu_tilde(:)-Mt_fort(:))))])

 
% 
% stop
% 

