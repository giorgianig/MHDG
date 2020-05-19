%-----------------------------------------------------
%  file : HDF5load.m
%  date : 01/04/2020
%  used to read HDF5 file
%  (developed by G.Giorgiani)
%-----------------------------------------------------
function res = HDF5load(hdf5_file)

info=h5info(hdf5_file);
% Datasets
for i=1:length(info.Datasets)
    name = info.Datasets(i).Name;
    data = getDataset('',name,hdf5_file);
    if nargout==0
        assignin('caller',name,data);
    else
        res.(name)=data;
    end
end
% Groups
for g = 1:length(info.Groups)
    [nameg,data] = getGroup(info.Groups(g),hdf5_file);
    if nargout==0
        assignin('caller',nameg,data);
    else
        res.(nameg)=data;
    end
    
end

function [nameg,res] = getGroup(info,hdf5_file)
p = strfind(info.Name,'/');
nameg =info.Name(p(end)+1:end);
for g = 1:length(info.Groups)
    [name,data] = getGroup(info.Groups(g),hdf5_file);
    res.(name) = data;
end
for i=1:length(info.Datasets)
    named = info.Datasets(i).Name;
    res.(named) = getDataset(info.Name,named,hdf5_file);
end

function data = getDataset(path,name,hdf5_file)
data=h5read(hdf5_file,[path,'/',name]);
