function M = readMatTxt(vec)

fid = fopen(vec,'r');
n1 = fscanf(fid,'%d',1);
n2 = fscanf(fid,'%d',1);
n3 = fscanf(fid,'%d',1);

if isempty(n2)
    n2 = 1;
end

if isempty(n3)
    n3=1;
end
nn = n1*n2*n3;
b = zeros(nn,1);

for i = 1:nn
    aux =fscanf(fid,'%f',1); 
    b(i) = aux;
end

fclose(fid);

M = reshape(b,[n1 n2 n3]);