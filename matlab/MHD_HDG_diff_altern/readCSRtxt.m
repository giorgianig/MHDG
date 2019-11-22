function Matrix = readCSRtxt(Mcsr)
% Read CSR matrix and convert to Matlab sparse

% open file
fid = fopen(Mcsr,'r');

% read n and nnz
n = fscanf(fid,'%d',1);
nnz = fscanf(fid,'%d',1);

col_ind = zeros(nnz,1);
val = col_ind;
row_ptr = zeros(n+1,1);
loc2glob = zeros(n,1);

% skip one
aux = fscanf(fid,'%s',1);
% read vals
for i = 1:nnz
    val(i) = fscanf(fid,'%f',1);
end

% skip one
aux = fscanf(fid,'%s',1);
% read cols
for i = 1:nnz
    col_ind(i) = fscanf(fid,'%d',1);
end

% skip one
aux = fscanf(fid,'%s',1);
% read rows
for i = 1:n+1
    row_ptr(i) = fscanf(fid,'%d',1);
end

% skip one
aux = fscanf(fid,'%s',1);
% read loc2glob
for i = 1:n
    loc2glob(i) = fscanf(fid,'%d',1);
end

% momentary fix sequential
% loc2glob = 1:n;

% close file
fclose(fid);

% Conversion
I = zeros(nnz,1);
J = I;
K = I;

for r = 1:numel(row_ptr)-1
    
    for i = row_ptr(r):(row_ptr(r+1)-1)
    
    I(i) = loc2glob(r);
    J(i) = col_ind(i);
    K(i) = val(i);
    end  
end

Matrix = sparse(I,J,K); 

