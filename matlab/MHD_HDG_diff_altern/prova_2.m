function B = prova_2(Bx,By)
n = size(Bx,1);
m = size(Bx,2);
B = zeros(size(Bx,1)*2,size(Bx,2)*4);
for j = 1:m
for i = 1:n

B((i-1)*2+1,(j-1)*4+1) = Bx(i,j);
B((i-1)*2+1,(j-1)*4+2) = By(i,j);
B((i-1)*2+2,(j-1)*4+3) = Bx(i,j);
B((i-1)*2+2,(j-1)*4+4) = By(i,j);

end
end