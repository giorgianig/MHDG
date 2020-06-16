function p = orthopoly3D(x,n)
% p = orthopoly3D(x,n)
% Computes the ortogonal base of 3D polynomials of degree less 
% or equal to n at the point x=(xi,eta,zeta) in the reference tetrahedra
%

xi = x(1); eta = x(2); zeta = x(3);

if (eta+zeta)==0 
    r = -1; s=1;
elseif zeta==1
    r = -1; s=1;  %or s=-1 (check that nothing changes)
else
    r = -2*(1+xi)/(eta+zeta)-1;
    s = 2*(1+eta)/(1-zeta)-1;
end
t = zeta;

p = orthopoly3D_rst([r,s,t],n);

