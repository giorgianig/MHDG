function [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
%
% [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

N = (n+1)*(n+2)/2 ;%number of nodes/polynomials
p = zeros(N,size(x,1));
dp_dxi  = zeros(N,size(x,1));
dp_deta = zeros(N,size(x,1));

r = x(:,1); s = x(:,2); 

xi = (1+r).*(1-s)/2-1;
eta = s;

dr_dxi  = 2./(1-eta);
dr_deta = 2*(1+xi)./(1-eta).^2;

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = ones(size(x,1),1); %1;
         q_i = ones(size(x,1),1); % 1;
         dp_i = zeros(size(x,1),1); %0;
         dq_i = zeros(size(x,1),1); %0;
     else
         p_i = jacobiPol(r,0,0,i);
         dp_i = jacobiPol(r,1,1,i-1)*(i+1)/2;    
         q_i = q_i.*(1-s)./2;
         dq_i = q_i.*(-i)./(1-s);
     end
     %Value for j
     j = nDeg-i;
     if j==0
        p_j = ones(size(x,1),1); %1;
        dp_j = zeros(size(x,1),1); %0; 
     else
        p_j = jacobiPol(s,2*i+1,0,j);
        dp_j = jacobiPol(s,2*i+2,1,j-1)*(j+2*i+2)/2;  
     end
     ncount= ncount+1;
     factor = sqrt( (2*i+1)*(i+j+1)/2 );
     %Normalized polinomial
     p(ncount,:)    = ( p_i.*q_i.*p_j )*factor;
     %Derivatives with respect to (r,s)
     dp_dr = ( (dp_i).*q_i.*p_j )*factor; 
     dp_ds = ( p_i.*(dq_i.*p_j+q_i.*dp_j) )*factor;
     %Derivatives with respect to (xi,eta)
     dp_dxi(ncount,:)  = dp_dr.*dr_dxi;
     dp_deta(ncount,:) = dp_dr.*dr_deta + dp_ds;
  end
end

