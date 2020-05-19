path2Fort = '/home/giorgio/Dropbox/Fortran/MHDG/test/';
M_fort = readMatTxt([path2Fort 'M.txt']);
B_fort = readMatTxt([path2Fort 'B.txt']);
C_fort = readMatTxt([path2Fort 'C.txt']);
Cdir_fort = readMatTxt([path2Fort 'Cdir.txt']);
Cv_fort = readMatTxt([path2Fort 'Cv.txt']);
H_fort = readMatTxt([path2Fort 'H.txt']);
Hdir_fort = readMatTxt([path2Fort 'Hdir.txt']);
D_fort = readMatTxt([path2Fort 'D.txt']);
E_fort = readMatTxt([path2Fort 'E.txt']);
P_fort = readMatTxt([path2Fort 'P.txt']);
Q_fort = readMatTxt([path2Fort 'Q.txt']);
Edir_fort = readMatTxt([path2Fort 'Edir.txt']);
S_fort = readMatTxt([path2Fort 'S.txt']);
Df_fort = readMatTxt([path2Fort 'Df.txt']);
Ef_fort = readMatTxt([path2Fort 'Ef.txt']);
Hf_fort = readMatTxt([path2Fort 'Hf.txt']);
Lf_fort = readMatTxt([path2Fort 'Lf.txt']);

Qf_fort = readMatTxt([path2Fort 'Qf.txt']);

fH_fort = readMatTxt([path2Fort 'fH.txt']);
U_fort = readMatTxt([path2Fort 'UU.txt']);
U0_fort = readMatTxt([path2Fort 'U0.txt']);


disp('M')
max(max(max(abs(M_fort-M))))
disp('B')
max(max(max(abs(B_fort-B))))
disp('C')
max(max(max(abs(C_fort-C))))
disp('Cdir')
max(max(max(abs(Cdir_fort-C_dir))))
disp('Cv')
max(max(max(abs(Cv_fort-Cv))))
disp('H')
max(max(max(abs(H_fort-H))))
disp('Hdir')
max(max(max(abs(Hdir_fort-Hdir))))
disp('D')
max(max(max(abs(D_fort-D-G))))
disp('E')
max(max(max(abs(E_fort-E))))

disp('P')
max(max(max(abs(P_fort-(P-Pb)))))
disp('Q')
max(max(max(abs(Q_fort-(Q-Qb)))))
disp('Lf')
max(max(max(abs(Lf_fort-(Lf)))))
disp('Qf')
max(max(max(abs(Qf_fort-(Qf)))))


disp('Edir')
max(max(max(abs(Edir_fort-Edir))))
disp('force')
max(max(max(abs(S_fort-force))))
disp('Df')
max(max(max(abs(Df_fort-Df))))
disp('Ef')
max(max(max(abs(Ef_fort-Ef))))
disp('Hf')
max(max(max(abs(Hf_fort-Hf))))
disp('fH')
max(max(max(abs(fH_fort-fH))))
disp('U')
max(max(max(abs(U_fort-UU))))
disp('U0')
max(max(max(abs(U0_fort-UU0))))