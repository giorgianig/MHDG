function C = expand_matrix(C11,C12,C21,C22)
n=size(C11,1);
m=size(C11,2);

C = zeros(2*size(C11,1),2*size(C11,2));
for i=1:n
    for j=1:m
        C((i-1)*2+1,(j-1)*2+1) = C((i-1)*2+1,(j-1)*2+1) + C11(i,j);
        C((i-1)*2+1,(j-1)*2+2) = C((i-1)*2+1,(j-1)*2+2) + C12(i,j);
        C((i-1)*2+2,(j-1)*2+1) = C((i-1)*2+2,(j-1)*2+1) + C21(i,j);
        C((i-1)*2+2,(j-1)*2+2) = C((i-1)*2+2,(j-1)*2+2) + C22(i,j);
    end
end
