clear
global  testcase axisym neq


testcase.n =50;
setneum   = 0; 
axisym      = 1;
neq           = 4;
dt =  7.2800e-1;
lscale = 1.901*1e-3;
% lscale =1;
nd              =2;
nf               = 3;
iel = 1902;
%% Load mesh & create reference element
loadData
Dirichlet_boundaries = setDirichletFaces(boundaryNames);
[F, F_dir, infoFaces, flipFace] = hdg_Preprocess(T,elementFaceInfo,boundaryNames,refEl,Dirichlet_boundaries);
Fcone = F(iel,:);

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
Np = size(refEl.NodesCoord,1);
np = size(refEl.NodesCoord1d,1);
perm = setperm(neq*np,neq);
for iface = 1:nf
    ind_v_L(iface,:) = neq*np*(iface-1)+ (1:neq*np);
end
ind_q  = (iel-1)*Np*neq*nd + (1:Np*neq*nd);
ind_u  = (iel-1)*Np*neq + (1:Np*neq);
ind_ut = reshape(bsxfun(@plus,(Fcone-1)*neq*np,(1:neq*np)'),neq*nf*np,1);

qe = q(ind_q);
ue = u(ind_u);
u0e = u0(ind_u);
ue_tilde = u_tilde(ind_ut);


flipFace_e = flipFace(iel,:);
for iface = 1:nf
    if flipFace_e(iface)
        ue_tilde(ind_v_L(iface,:)) = ue_tilde(ind_v_L(iface,perm));
    end
end
ind = 3:neq:(Np*neq-1);

sum(inv(iL)*qe+B*ue-C*ue_tilde)

res = (M*(ue-u0e)/dt - Cv*ue+H*ue_tilde+D*ue-E*ue_tilde+TQ*qe+P*qe-(S+Tf))
sum(res)


enInc = lscale^3*M(3:neq:end-1,:)*(ue-u0e)
% 
% M(ind_ul,ind_ul)*(u(ind_u,3) -u0(ind_u,3))...
% -Cv(ind_ul,ind_ul)*u(ind_u,3) + H(ind_ul,ind_utl)*u_tilde(ind_ut,3)+D(ind_ul,ind_ul)*u(ind_u,3) -E(ind_ul,ind_utl)*u_tilde(ind_ut,3)+...
%     TQ(ind_ul,ind_ql)*col(transpose(q(ind_q,5:6)))+P(ind_ul,ind_ql)*col(transpose(q(ind_q,5:6)))-(S(ind_ul) +Tf(ind_ul))


aux = - (- Cv*ue+H*ue_tilde+D*ue-E*ue_tilde+TQ*qe+P*qe-(S+Tf));

aux_1 = -(- Cv*ue+H*ue_tilde + TQ*qe-Tf);

aux_2 = -(D*ue-E*ue_tilde);

aux_3 = -P*qe;

disp(['sum(aux(ind))= ', num2str(sum(aux(ind)))])
disp(['sum(aux_1(ind))= ', num2str(sum(aux_1(ind)))])
disp(['sum(aux_2(ind))= ', num2str(sum(aux_2(ind)))])
disp(['sum(aux_3(ind))= ', num2str(sum(aux_3(ind)))])



