function M = readMatTxt(vec)

fid = fopen(vec,'r');
n1 = fscanf(fid,'%d',1);
n2 = fscanf(fid,'%d',1);
n3 = fscanf(fid,'%d',1);
nn = n1*n2*n3;
b = zeros(nn,1);

for i = 1:nn
    b(i) = fscanf(fid,'%f',1);
end

fclose(fid);

M = reshape(b,[n1 n2 n3]);