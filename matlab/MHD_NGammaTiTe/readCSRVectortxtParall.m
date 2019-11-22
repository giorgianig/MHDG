function Vec = readCSRVectortxtParall(Mcsr,nproc)
% Read CSR matrix and convert to Matlab sparse

dimension = 0;
Loc2glob = [];
Vals     = [];
for iproc = 1:nproc
    % open file
    fid = fopen([Mcsr '_' num2str(iproc-1) '_' num2str(nproc) '.txt'],'r');
    
    % read n and nnz
    n = fscanf(fid,'%d',1);
    
    val = zeros(n,1);
    loc2glob = zeros(n,1);
    
    % skip one
    aux = fscanf(fid,'%s',1);
    % read vals
    for i = 1:n
        val(i) = fscanf(fid,'%f',1);
    end
       
    % skip one
    aux = fscanf(fid,'%s',1);
    % read loc2glob
    for i = 1:n
        loc2glob(i) = fscanf(fid,'%d',1);
    end
      
    % close file
    fclose(fid);
    
    dimension = dimension+n;   
    
    Loc2glob = [Loc2glob;loc2glob];
    Vals     = [Vals;val];
end
Vec = zeros(dimension,1);
Vec(Loc2glob) = Vals;
