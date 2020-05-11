function Vec = readCSRVectortxt(Mcsr)
% Read CSR matrix and convert to Matlab sparse

% open file
fid = fopen(Mcsr,'r');

% read n and nnz
n = fscanf(fid,'%d',1);

val = zeros(n,1);

% skip one
aux = fscanf(fid,'%s',1);

% read vals
for i = 1:n
    val(i) = fscanf(fid,'%f',1);
end

% close file
fclose(fid);

Vec     = val;
