%% Debug fortran
% M_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/M.txt');
% Cv_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Cv.txt');
% H_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/H.txt');
% if any(Dirichlet_boundaries)
%     Hdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Hdir.txt');
%     Edir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Edir.txt');
%     Thdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Thdir.txt');
%     Cdir_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Cdir.txt');
% end
% D_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/D.txt');
% E_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/E.txt');
% S_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/S.txt');
% UU_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/UU.txt');
% U0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/U0.txt');
% Hf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Hf.txt');
% Df_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Df.txt');
% Ef_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Ef.txt');
% fH_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/fH.txt');
% B_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/B.txt');
% C_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/C.txt');
% P_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/P.txt');
% G_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/G.txt');
% iL_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/iL.txt');
% Lf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Lf.txt');
% Qf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Qf.txt');
% LL_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/LL.txt');
% L0_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/L0.txt');
% TQ_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/TQ.txt');
% TQhf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/TQhf.txt');
% Tfhf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Tfhf.txt');
% Tf_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Tf.txt');
% matK= readCSRtxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/MatK.txt');
% f_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/f.txt');




M_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/M.h5');
Cv_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Cv.h5');
H_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/H.h5');
if any(Dirichlet_boundaries)
    Hdir_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Hdir.h5');
    Edir_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Edir.h5');
    Thdir_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Thdir.h5');
    Cdir_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Cdir.h5');
end
D_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/D.h5');
E_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/E.h5');
S_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/S.h5');
UU_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/UU.h5');
U0_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/U0.h5');
Hf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Hf.h5');
Df_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Df.h5');
Ef_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Ef.h5');
fH_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/fH.h5');
B_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/B.h5');
C_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/C.h5');
P_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/P.h5');
G_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/G.h5');
iL_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/iL.h5');
Lf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Lf.h5');
Qf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Qf.h5');
LL_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/LL.h5');
L0_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/L0.h5');
TQ_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/TQ.h5');
TQhf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/TQhf.h5');
Tfhf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Tfhf.h5');
Tf_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/Tf.h5');
matK= HDF5readCSRmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/matK.h5');
f_fort = HDF5readmat('/home/giorgio/Dropbox/Fortran/MHDG_3D/test/rhs.h5');




if (max(max(max(abs(M_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(M_fort)))));
end
err = max(max(max(abs(M-M_fort))))/ref;
disp(['Error in M: ', num2str(err)])




if (max(max(max(abs(Cv_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Cv_fort)))));
end
err = max(max(max(abs(Cv-TU-Cv_fort))))/ref;
disp(['Error in Cv: ', num2str(err)])





if (max(max(max(abs(H_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(H_fort)))));
end
err = max(max(max(abs(H-TUh-H_fort))))/ref;
disp(['Error in H: ', num2str(err)])





if (max(max(max(abs(Hf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Hf_fort)))));
end
err = max(max(max(abs(Hf-TUhf-Hf_fort))))/ref;
disp(['Error in Hf: ', num2str(err)])






if any(Dirichlet_boundaries)    
    if (max(max(max(abs(Hdir_fort)))))==0
        ref = 1;
    else
        ref = (max(max(max(abs(Hdir_fort)))));
    end
    err = max(max(max(abs(Hdir-Hdir_fort))))/ref;
    disp(['Error in Hdir: ', num2str(err)])
end






if (max(max(max(abs(D_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(D_fort)))));
end
err = max(max(max(abs(D-D_fort))))/ref;
disp(['Error in D: ', num2str(err)])







if (max(max(max(abs(E_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(E_fort)))));
end
err = max(max(max(abs(E-E_fort))))/ref;
disp(['Error in E: ', num2str(err)])







if any(Dirichlet_boundaries)

if (max(max(max(abs(Edir_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Edir_fort)))));
end    
    err = max(max(max(abs(Edir-Edir_fort))))/ref;
    disp(['Error in Edir: ', num2str(err)])
end







if (max(max(max(abs(S_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(S_fort)))));
end
err = max(max(max(abs(force-S_fort))))/ref;
disp(['Error in f: ', num2str(err)])







if (max(max(max(abs(Df_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Df_fort)))));
end
err = max(max(max(abs(Df-Df_fort))))/ref;
disp(['Error in Df: ', num2str(err)])






if (max(max(max(abs(Ef_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Ef_fort)))));
end
err = max(max(max(abs(Ef-Ef_fort))))/ref;
disp(['Error in Ef: ', num2str(err)])







if (max(max(max(abs(fH_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(fH_fort)))));
end
err = max(max(max(abs(fH-fH_fort))))/ref;
disp(['Error in fH: ', num2str(err)])







if (max(max(max(abs(B_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(B_fort)))));
end
err = max(max(max(abs(B-B_fort))))/ref;
disp(['Error in B: ', num2str(err)])








if (max(max(max(abs(C_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(C_fort)))));
end
err = max(max(max(abs(C-C_fort))))/ref;
disp(['Error in C: ', num2str(err)])









if any(Dirichlet_boundaries)
if (max(max(max(abs(Cdir_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Cdir_fort)))));
end    
    err = max(max(max(abs(C_dir-Cdir_fort))))/ref;
    disp(['Error in Cdir: ', num2str(err)])
end







if (max(max(max(abs(P_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(P_fort)))));
end
err = max(max(max(abs(P-Pb-Q+Qb-P_fort))))/ref;
disp(['Error in P: ', num2str(err)])

% err = max(max(max(abs(G-G_fort))))/max(max(max(abs(G_fort))));
% disp(['Error in G: ', num2str(err)])








if (max(max(max(abs(iL_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(iL_fort)))));
end
err = max(max(max(abs(invL-iL_fort))))/ref;
disp(['Error in iL: ', num2str(err)])








if (max(max(max(abs(Lf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Lf_fort)))));
end
err = max(max(max(abs(Lf-Lf_fort))))/ref;
disp(['Error in Lf: ', num2str(err)])








if (max(max(max(abs(Qf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Qf_fort)))));
end
err = max(max(max(abs(Qf-Qf_fort))))/ref;
disp(['Error in Qf: ', num2str(err)])







if (max(max(max(abs(TQ_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(TQ_fort)))));
end
err = max(max(max(abs(TQ-TQh-TQ_fort))))/ref;
disp(['Error in TQ: ', num2str(err)])







if (max(max(max(abs(TQhf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(TQhf_fort)))));
end
err = max(max(max(abs(TQhf-TQhf_fort))))/ref;
disp(['Error in TQf: ', num2str(err)])







if (max(max(max(abs(Tfhf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Tfhf_fort)))));
end
err = max(max(max(abs(Tfhf-Tfhf_fort))))/ref;
disp(['Error in Tfhf: ', num2str(err)])







if (max(max(max(abs(Tf_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Tf_fort)))));
end
err = max(max(max(abs(Tf-Tfh-Tf_fort))))/ref;
disp(['Error in Tf: ', num2str(err)])









if any(Dirichlet_boundaries)
if (max(max(max(abs(Thdir_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(Thdir_fort)))));
end    
    err = max(max(max(abs(TUhdir-Thdir_fort))))/ref;
    disp(['Error in Thdir: ', num2str(err)])
end







if (max(max(max(abs(UU_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(UU_fort)))));
end
err = max(max(max(abs(UU-UU_fort))))/ref;
disp(['Error in UU: ', num2str(err)])










if (max(max(max(abs(U0_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(U0_fort)))));
end
err = max(max(max(abs(UU0-U0_fort))))/ref;
disp(['Error in U0: ', num2str(err)])








if (max(max(max(abs(LL_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(LL_fort)))));
end
err = max(max(max(abs(LL-LL_fort))))/ref;
disp(['Error in LL: ', num2str(err)])








if (max(max(max(abs(L0_fort)))))==0
    ref = 1;
else
    ref = (max(max(max(abs(L0_fort)))));
end
err = max(max(max(abs(LL0-L0_fort))))/ref;
disp(['Error in L0: ', num2str(err)])








disp(['Error in global matrix: ' num2str(max(max(abs(matK-K))))])
disp(['Error in global RHS: ' num2str((max(abs(f-f_fort))))])


















