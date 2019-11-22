function p = orthopoly2D(x,n)
% p = orthopoly2D(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(xi,eta) in the reference triangle
%

xi = x(:,1); eta = x(:,2); 

r = 2*(1+xi)./(1-eta)-1;
s = eta;
arrange=find(eta==1);
r(arrange)=-1;
s(arrange)=1;


% if (r~=x(1) || s~=x(2)), stop,end

p = orthopoly2D_rst([r,s],n);