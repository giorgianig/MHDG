home
clear
close all
Mat_par = readCSRtxtParall('/home/giorgio/Dropbox/Fortran/MHDG_3D_ptor/test/MatK',8);
Mat_seq = readCSRtxt('/home/giorgio/Dropbox/Fortran/MHDG_3D_ptor/test/MatK_seq.txt');
size(Mat_par)
size(Mat_seq)
max(max(abs(Mat_seq-Mat_par)))
rhs_seq = readCSRVectortxt('/home/giorgio/Dropbox/Fortran/MHDG_3D_ptor/test/rhs_seq.txt');
rhs_par = readCSRVectortxtParall('/home/giorgio/Dropbox/Fortran/MHDG_3D_ptor/test/rhs',8);
max(max(abs(rhs_par-rhs_seq)))
% sol_ser = Mat_seq\rhs_seq;
% sol_par = Mat_par\rhs_par;
% max(max(abs(sol_ser-sol_par)));
% max(max(abs(sol_ser-sol_par)))