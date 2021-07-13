function res=HDF5readmat(hdf5_file)

info=hdf5info(hdf5_file);
name=info.GroupHierarchy.Datasets.Name(2:end);
data=hdf5read(hdf5_file,name);
res=data;
% for i=1:length(info.GroupHierarchy.Datasets)
%  name{i}=info.GroupHierarchy.Datasets(i).Name(2:end);
%  if (length(info.GroupHierarchy.Datasets(i).Dims)==1)
%    data{i}=hdf5read(hdf5_file,name{i})';
%  else
%    data{i}=hdf5read(hdf5_file,name{i});
%  end
%  assignin('caller',name{i},data{i}) 
% end