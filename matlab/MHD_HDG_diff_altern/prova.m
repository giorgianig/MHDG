function B = prova(Bx,By)
n = size(Bx,1);
m = size(Bx,2);
B = zeros(size(Bx,1)*4,size(Bx,2)*2);

for j = 1:m
for i = 1:n

B((i-1)*4+1,(j-1)*2+1) = Bx(i,j);
B((i-1)*4+2,(j-1)*2+1) = By(i,j);
B((i-1)*4+3,(j-1)*2+2) = Bx(i,j);
B((i-1)*4+4,(j-1)*2+2) = By(i,j);

end
end