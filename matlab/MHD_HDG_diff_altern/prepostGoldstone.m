% pre-postprocess for Goldston

path2fort = '/home/giorgio/Goldston/YesDrift/M2';
solname = 'Sol_Circle_LIM_InfThin_SmallSOL0.77_refCorn0.001_refSol0.008_DoubleSol_P5_Diff.25000E-03.h5';

% start
pos = strfind(solname,'_P');
for i=1:10
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

save([path2save, solname(1:end-3),'.mat'],'u','u_tilde')