%% Debug fortran
disp('Reading fortran matrices')
path2fort = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/';
M_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/M.txt');
Cv_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Cv.txt');
H_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/H.txt');
if any(Dirichlet_boundaries)
    Hdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Hdir.txt');
    Edir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Edir.txt');
%     Thdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Thdir.txt');
    Cdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Cdir.txt');
end
D_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/D.txt');
E_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/E.txt');
S_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/S.txt');
UU_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/UU.txt');
U0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/U0.txt');
Hf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Hf.txt');
Df_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Df.txt');
Ef_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Ef.txt');
fH_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/fH.txt');
B_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/B.txt');
C_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/C.txt');
P_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/P.txt');
G_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/G.txt');
iL_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/iL.txt');
Lf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Lf.txt');
Qf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Qf.txt');
LL_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/LL.txt');
L0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/L0.txt');
% TQ_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/TQ.txt');
% TQhf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/TQhf.txt');
% Tfhf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Tfhf.txt');
% Tf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/Tf.txt');
% matK= readCSRtxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/MatK.txt');
% f_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test_anal/f.txt');


err = max(max(max(abs(M-M_fort))));
disp(['Error in M: ', num2str(err)])

err = max(max(max(abs(Cv-Cv_fort))));
disp(['Error in Cv: ', num2str(err)])

err = max(max(max(abs(H-H_fort))));
disp(['Error in H: ', num2str(err)])

err = max(max(max(abs(Hf-Hf_fort))));
disp(['Error in Hf: ', num2str(err)])

if any(Dirichlet_boundaries)
    err = max(max(max(abs(Hdir-Hdir_fort))));
    disp(['Error in Hdir: ', num2str(err)])
end

err = max(max(max(abs(D-D_fort))));
disp(['Error in D: ', num2str(err)])

err = max(max(max(abs(E-E_fort))));
disp(['Error in E: ', num2str(err)])

if any(Dirichlet_boundaries)
    err = max(max(max(abs(Edir-Edir_fort))));
    disp(['Error in Edir: ', num2str(err)])
end

err = max(max(max(abs(force-S_fort))));
disp(['Error in f: ', num2str(err)])

err = max(max(max(abs(Df-Df_fort))));
disp(['Error in Df: ', num2str(err)])

err = max(max(max(abs(Ef-Ef_fort))));
disp(['Error in Ef: ', num2str(err)])

err = max(max(max(abs(fH-fH_fort))));
disp(['Error in fH: ', num2str(err)])

err = max(max(max(abs(B-B_fort))));
disp(['Error in B: ', num2str(err)])

err = max(max(max(abs(C-C_fort))));
disp(['Error in C: ', num2str(err)])

if any(Dirichlet_boundaries)
    err = max(max(max(abs(C_dir-Cdir_fort))));
    disp(['Error in Cdir: ', num2str(err)])
end

err = max(max(max(abs(P-Pb-Q+Qb-P_fort))));
disp(['Error in P: ', num2str(err)])

err = max(max(max(abs(G-G_fort))));
disp(['Error in G: ', num2str(err)])

err = max(max(max(abs(invL-iL_fort))));
disp(['Error in iL: ', num2str(err)])

err = max(max(max(abs(Lf-Lf_fort))));
disp(['Error in Lf: ', num2str(err)])

err = max(max(max(abs(Qf-Qf_fort))));
disp(['Error in Qf: ', num2str(err)])

% err = max(max(max(abs(TQ-TQh-TQ_fort))));
% disp(['Error in TQ: ', num2str(err)])
% 
% err = max(max(max(abs(TQhf-TQhf_fort))));
% disp(['Error in TQf: ', num2str(err)])
% 
% err = max(max(max(abs(Tfhf-Tfhf_fort))));
% disp(['Error in Tfhf: ', num2str(err)])
% 
% err = max(max(max(abs(Tf-Tfh-Tf_fort))));
% disp(['Error in Tf: ', num2str(err)])
% 
% if any(Dirichlet_boundaries)
%     err = max(max(max(abs(TUhdir-Thdir_fort))));
%     disp(['Error in Thdir: ', num2str(err)])
% end
err = max(max(max(abs(UU-UU_fort))));
disp(['Error in UU: ', num2str(err)])

err = max(max(max(abs(UU0-U0_fort))));
disp(['Error in U0: ', num2str(err)])

err = max(max(max(abs(LL-LL_fort))));
disp(['Error in LL: ', num2str(err)])

err = max(max(max(abs(LL0-L0_fort))));
disp(['Error in L0: ', num2str(err)])

% disp(['Error in global matrix: ' num2str(max(max(abs(matK-K))))])
% disp(['Error in global RHS: ' num2str((max(abs(f-f_fort))))])


















