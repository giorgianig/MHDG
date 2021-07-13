%% Debug fortran
clear

%**********************************
% Parallel/serial
%**********************************
nproc=16; % number of MPI task for parallel (mesh separation)
solpath = '/home/bluce/MHDG_sim/West1090_P4/2021_07_01_NGT2D/';
save = 1; % True to save only I,J,K then use python script to compute on a cluster the global matrix (useful for 3D matrix for whom condest crashes)

HDF5readCSRmat([solpath,'Mat'],nproc,save)


% I = [];
% J = [];
% K = [];
% 
% for iproc = 1:nproc
%     if nproc == 1
%         hdf5_file = [solpath,'Mat','_ijk.h5'];
%     else
%         hdf5_file = [solpath,'Mat','_ijk_',num2str(iproc),'_',num2str(nproc),'.h5'];
%     end
%     
%     info=hdf5info(hdf5_file);
%     Iloc = hdf5read(hdf5_file,'I');
%     Jloc = hdf5read(hdf5_file,'J');
%     Kloc = hdf5read(hdf5_file,'K');
%     
%     I = [I; Iloc];
%     J = [J; Jloc];
%     K = [K; Kloc];
%     
% end
% 
% Matrix = sparse(I(I~=0),J(I~=0),K(I~=0));
% disp(['Error in global matrix: ', num2str(condest(Matrix,1),'%e')])

