function HDF5readCSRmat(hdf5_file_gen,nproc,save)
if nargin == 2
    save = 0;
end

dimension = 0;
I = [];
J = [];
K = [];

for iproc = 1:nproc
    if nproc == 1
        hdf5_file = [hdf5_file_gen,'.h5'];
        hdf5_save = [hdf5_file_gen,'_ijk.h5'];
    else
        hdf5_file = [hdf5_file_gen,'_',num2str(iproc),'_',num2str(nproc),'.h5'];
        hdf5_save = [hdf5_file_gen,'_ijk_',num2str(iproc),'_',num2str(nproc),'.h5'];
    end
    
    info=hdf5info(hdf5_file);
    MatCSR.n = hdf5read(hdf5_file,'n');
    MatCSR.nnz = hdf5read(hdf5_file,'nnz');
    MatCSR.cols = hdf5read(hdf5_file,'cols');
    MatCSR.rowptr = hdf5read(hdf5_file,'rowptr');
    MatCSR.loc2glob = hdf5read(hdf5_file,'loc2glob');
    MatCSR.vals = hdf5read(hdf5_file,'vals');
    
    Iloc = zeros(MatCSR.nnz,1);
    Jloc = Iloc;
    Kloc = Iloc;
    for r = 1:numel(MatCSR.rowptr)-1     
        for i = MatCSR.rowptr(r):(MatCSR.rowptr(r+1)-1)
            Iloc(i) = MatCSR.loc2glob(r);
            Jloc(i) = MatCSR.cols(i);
            Kloc(i) = MatCSR.vals(i);
        end
    end
    if save
        h5create(hdf5_save,'/I',size(Iloc))
        h5write(hdf5_save,'/I',Iloc)
        h5create(hdf5_save,'/J',size(Jloc))
        h5write(hdf5_save,'/J',Jloc)
        h5create(hdf5_save,'/K',size(Kloc))
        h5write(hdf5_save,'/K',Kloc)
    else
        I = [I; Iloc];
        J = [J; Jloc];
        K = [K; Kloc];
        dimension = dimension+MatCSR.n;
    end
end
if ~save
    Matrix = sparse(I(I~=0),J(I~=0),K(I~=0));
    disp(['Condition number: ', num2str(condest(Matrix),'%e')])
end

end
