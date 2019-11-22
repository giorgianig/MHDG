clear
clc

Cv_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Cv_seq.txt');
Cv_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Cv_par.txt');
res = max(abs(Cv_seq(:)-Cv_par(:)));
disp(['Error in Cv: ', num2str(res)])
if res<1e-5,clear Cv_seq Cv_par,end


Hf_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Hf_seq.txt');
Hf_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Hf_par.txt');
res =max(abs(Hf_seq(:)-Hf_par(:)));
disp(['Error in Hf: ', num2str(res)])
if res<1e-5,clear Hf_seq Hf_par,end

H_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/H_seq.txt');
H_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/H_par.txt');
res = max(abs(H_seq(:)-H_par(:)));
disp(['Error in H: ', num2str(res)])
if res<1e-5,clear H_seq H_par,end

E_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/E_seq.txt');
E_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/E_par.txt');
res = max(abs(E_seq(:)-E_par(:)));
disp(['Error in E: ', num2str(res)])
if res<1e-5,clear E_seq E_par,end


D_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/D_seq.txt');
D_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/D_par.txt');
res = max(abs(D_seq(:)-D_par(:)));
disp(['Error in D: ', num2str(res)])
if res<1e-5,clear D_seq D_par,end


Ef_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Ef_seq.txt');
Ef_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Ef_par.txt');
res = max(abs(Ef_seq(:)-Ef_par(:)));
disp(['Error in Ef: ', num2str(res)])
if res<1e-5,clear Ef_seq Ef_par,end


Df_seq = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Df_seq.txt');
Df_par = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/Df_par.txt');
res = max(abs(Df_seq(:)-Df_par(:)));
disp(['Error in Df: ', num2str(res)])
if res<1e-5,clear Df_seq Df_par,end




