function [p,dp_dxi,dp_deta,dp_dzet] = orthopoly3D_deriv_rst(x,n)
%
% [p,dp_dxi,dp_deta,dp_dzet] = orthopoly3D_deriv_rst(x,n)
% Computes the ortogonal base of 3D polynomials of degree less 
% or equal to n at the point x=(r,s,t) in [-1,1]^3
%

N = (n+1)*(n+2)*(n+3)/6 ;%number of nodes/polynomials
p = zeros(N,1);
dp_dxi  = zeros(N,1);
dp_deta = zeros(N,1);
dp_dzet = zeros(N,1);

r = x(1); s = x(2); t = x(3);

eta = (1/2)*(s-s*t-1-t);
xi = -(1/2)*(r+1)*(eta+t)-1;
zeta = t;

dr_dxi  = -2/(eta+zeta);
dr_deta = 2*(1+xi)/(eta+zeta)^2;
dr_dzet = dr_deta;

ds_deta = 2/(1-zeta);
ds_dzet = 2*(1+eta)/(1-zeta)^2;

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = 1;  q_i = 1;  dp_i = 0; dq_i = 0;
     else
         p_i = jacobiPol(r,0,0,i); dp_i = jacobiPol(r,1,1,i-1)*(i+1)/2;    
         q_i = q_i*(1-s)/2; dq_i = q_i*(-i)/(1-s);
     end
     %Loop increasing j
     for j = 0:(nDeg-i)
        if j==0
           p_j = 1;  q_j = ((1-t)/2)^i;  dp_j = 0; dq_j = q_j*(-(i+j))/(1-t);  
        else
           p_j = jacobiPol(s,2*i+1,0,j); dp_j = jacobiPol(s,2*i+2,1,j-1)*(j+2*i+2)/2;  
           q_j = q_j*(1-t)/2;  dq_j = q_j*(-(i+j))/(1-t);
        end
        %Value for k
        k = nDeg-(i+j);
        if k==0
           p_k = 1;  dp_k = 0; 
        else
           p_k = jacobiPol(t,2*(i+j)+2,0,k);  dp_k = jacobiPol(t,2*(i+j)+3,1,k-1)*(k+2*i+2*j+3)/2;
        end
        ncount= ncount+1;
        factor = sqrt( (2*i+1)*(i+j+1)*(2*(i+j+k)+3)/4 );
        %Normalized polinomial
        p(ncount)    = ( p_i*q_i*p_j*q_j*p_k )*factor;
        %Derivatives with respect to (r,s,t)
        dp_dr = ( (dp_i)*q_i*p_j*q_j*p_k )*factor; 
        dp_ds = ( p_i*(dq_i*p_j+q_i*dp_j)*q_j*p_k )*factor;
        dp_dt = ( p_i*q_i*p_j*(dq_j*p_k+q_j*dp_k) )*factor;
        %Derivatives with respect to (xi,eta,zeta)
        dp_dxi(ncount)  = dp_dr*dr_dxi;
        dp_deta(ncount) = dp_dr*dr_deta + dp_ds*ds_deta;
        dp_dzet(ncount) = dp_dr*dr_dzet + dp_ds*ds_dzet + dp_dt;
    end
  end
end