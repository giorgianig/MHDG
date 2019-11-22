function res=HDF5readCSRmat(hdf5_file)

info=hdf5info(hdf5_file);
MatCSR.n = hdf5read(hdf5_file,'n');
MatCSR.nnz = hdf5read(hdf5_file,'nnz');
MatCSR.cols = hdf5read(hdf5_file,'cols');
MatCSR.rowptr = hdf5read(hdf5_file,'rowptr');
MatCSR.loc2glob = hdf5read(hdf5_file,'loc2glob');
MatCSR.vals = hdf5read(hdf5_file,'vals');

res = convertMatrixCSRtoIJK(MatCSR);