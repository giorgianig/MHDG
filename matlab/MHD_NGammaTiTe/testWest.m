clear
global  testcase axisym neq


testcase.n =2;
setneum   = 0; 
axisym      = 1;
neq           = 4;
dt =  7.2800e-1;
lscale = 1.901*1e-3;
lscale =1;
nd              =2;
nf               = 3;
iel = 350;
%% Load mesh & create reference element
loadData
Dirichlet_boundaries = setDirichletFaces(boundaryNames);
[F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,refEl,Dirichlet_boundaries);

%% Reading Fortran Matrices
M    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/M.txt');
Cv   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Cv.txt');
H     = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/H.txt');
D     = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/D.txt');
E     = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/E.txt');
S     = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/S.txt');
UU  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/UU.txt');
U0  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/U0.txt');
Hf   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Hf.txt');
Df   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Df.txt');
Ef    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Ef.txt');
fH   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/fH.txt');
B    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/B.txt');
C    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/C.txt');
P    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/P.txt');
G   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/G.txt');
iL   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/iL.txt');
Lf   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Lf.txt');
Qf  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Qf.txt');
LL  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/LL.txt');
L0  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/L0.txt');
TQ  = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/TQ.txt');
TQhf = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/TQhf.txt');
Tfhf   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Tfhf.txt');
Tf   = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Tf.txt');
q    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/q.txt');
u    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/u.txt');
u0    = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/u0.txt');
u_tilde = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/u_tilde.txt');

%% 
Nu  = numel(T);
Nut = numel(u_tilde);
q  = transpose(reshape(q,[neq*nd,Nu]));
u = transpose(reshape(u,[neq,Nu]));
u0 = transpose(reshape(u0,[neq,Nu]));
u_tilde = transpose(reshape(u_tilde,[neq,Nut/neq]));
Np = size(refEl.NodesCoord,1);
np = size(refEl.NodesCoord1d,1);
ind_q  = (iel-1)*Np + (1:Np);
ind_u  = (iel-1)*Np + (1:Np);
ind_ut = bsxfun(@plus,(F(iel,:)-1)*np,(1:np)');
for ifa =1:nf
    if flipFace(iel,ifa)
    ind_ut(:,ifa) = flipud(ind_ut(:,ifa));
    end
end
ind_ut = col(ind_ut);

ind_ql = col([5:neq*nd:(neq*Np*nd-3);6:neq*nd:(neq*Np*nd-2)]);
ind_ul  = 3:neq:(neq*Np-1);
ind_utl = 3:neq:(neq*np*nf-1);



u(ind_u,3)

eninc = lscale^3*sum(M(ind_ul,ind_ul)*(u(ind_u,3) -u0(ind_u,3)))

intsol = computeIntSol(X*lscale,T,u -u0,refEl);
intsol(iel,3)





M(ind_ul,ind_ul)*(u(ind_u,3) -u0(ind_u,3))...
-Cv(ind_ul,ind_ul)*u(ind_u,3) + H(ind_ul,ind_utl)*u_tilde(ind_ut,3)+D(ind_ul,ind_ul)*u(ind_u,3) -E(ind_ul,ind_utl)*u_tilde(ind_ut,3)+...
    TQ(ind_ul,ind_ql)*col(transpose(q(ind_q,5:6)))+P(ind_ul,ind_ql)*col(transpose(q(ind_q,5:6)))-(S(ind_ul) +Tf(ind_ul))




