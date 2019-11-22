% plot Fortran solution
clear

solpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/';
solname = 'Sol_Circle_limInfThin_h0.1_P8_Diff.10000E+02.h5';
path2save =    '/home/giorgio/Dropbox/Matlab/MHD/MHD_NGammaEnergy/Saves/';

% start
pos = strfind(solname,'_P');
for i=1:10;
    if strcmp(solname(pos+i),'_')
        pos = pos+i-1;
        break
    end
end
meshname = [solname(5:pos) '.h5'];

HDF5load([solpath,meshname])
HDF5load([solpath,solname])
u = u';
u_tilde = u_tilde';
gradU = q';

save([path2save, solname(1:end-3),'.mat'],'u','u_tilde','gradU')