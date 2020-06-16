function [p,dp] = orthopoly1D_deriv(x,n)
% [p,dp] = orthopoly1D_deriv(x,n)
% Computes the ortogonal base of 1D polynomials of degree less 
% or equal to n at the point x in [-1,1]
%

p = zeros(n+1,size(x,1));
dp = zeros(n+1,size(x,1));

p(1,:) = 1/sqrt(2);
dp(1,:) = 0;
for i=1:n
    factor = sqrt((2*i+1)/2);
    p(i+1,:)  = jacobiPol(x,0,0,i)*factor; 
    dp(i+1,:) = jacobiPol(x,1,1,i-1)*((i+1)/2)*factor;
end
