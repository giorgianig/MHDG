function p = orthopoly1D(x,n)
% p = orthopoly1D(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x in [-1,1]
%

p = zeros(n+1,size(x,1));

for i=0:n
    p(i+1,:) = jacobiPol(x,0,0,i)*sqrt((2*i+1)/2); 
end
