% Convert a mesh in .mat to .hdf5 

fileName = [fileName(1:end-4),'.h5'];


if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end
    
% Starting conversion
delete(fileName)
Nnodes = size(X,1);
Nelems = size(T,1);
Ndim   = size(X,2);
Nnodesperelem = size(T,2);
h5create(fileName,'/T',size(T));
h5create(fileName,'/X',size(X));
h5create(fileName,'/Nnodes',1);
h5create(fileName,'/Ndim',1);
h5create(fileName,'/Nelems',1);
h5create(fileName,'/elemType',1);
h5create(fileName,'/Nnodesperelem',1);
% Boundaries
aux = whos('Tb*');
boundaryNames = {aux.name};
clear aux
nOfBound = numel(boundaryNames);
s1 = 0;
boundInt = zeros(nOfBound,1);
for ibound = 1:nOfBound
    name = boundaryNames{ibound};
    eval(['aux_bound = Tb_' name(4:end) ';'])
    s1 = s1+size(aux_bound,1);
    if any(strcmpi(name,{'Tb_Dirichlet','Tb_Diriclet'}))
        boundInt(ibound) = 1;
    elseif strcmpi(name,'Tb_LEFT')
        boundInt(ibound) = 2;
    elseif strcmpi(name,'Tb_RIGHT')
        boundInt(ibound) = 3;
    elseif strcmpi(name,'Tb_UP')
        boundInt(ibound) = 4;
    elseif strcmpi(name,'Tb_DOWN')
        boundInt(ibound) = 5;
    elseif strcmpi(name,'Tb_WALL')
        boundInt(ibound) = 6;
    elseif strcmpi(name,'Tb_LIM')
        boundInt(ibound) = 7; 
    elseif strcmpi(name,'Tb_IN')
        boundInt(ibound) = 8; 
    elseif strcmpi(name,'Tb_OUT')
        boundInt(ibound) = 9;
    elseif strcmpi(name,'Tb_ULIM')        
        boundInt(ibound) = 10;
    end
end
s2 = size(aux_bound,2);
Tb = zeros(s1,s2);
boundaryFlag = zeros(s1,1);
s1 = 1;
for ibound = 1:nOfBound
    name = boundaryNames{ibound};
    eval(['aux_bound = Tb_' name(4:end) ';'])
    Tb(s1:size(aux_bound,1)+s1-1,:) = aux_bound;
    boundaryFlag(s1:size(aux_bound,1)+s1-1) = boundInt(ibound);
    s1=s1+size(aux_bound,1);
end
h5create(fileName,'/Nnodesperface',1);
h5create(fileName,'/Nextfaces',1);
h5create(fileName,'/Tb',size(Tb));
h5create(fileName,'/boundaryFlag',size(boundaryFlag));

% Write in the database
h5disp(fileName);
h5write(fileName,'/T',T);
h5write(fileName,'/X',X);
h5write(fileName,'/Nnodes',Nnodes);
h5write(fileName,'/Ndim',Ndim);
h5write(fileName,'/Nelems',Nelems);
h5write(fileName,'/Nnodesperelem',Nnodesperelem);
h5write(fileName,'/elemType',elemType);
h5write(fileName,'/Nextfaces',size(Tb,1));
h5write(fileName,'/Nnodesperface',s2);
h5write(fileName,'/Tb',Tb);
h5write(fileName,'/boundaryFlag',boundaryFlag);

% Check
h5disp(fileName);
% h5read(fileName,'/Ndim')