function p = orthopoly2D_rst(x,n)
% p = orthopoly2D_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

N = (n+1)*(n+2)/2 ;%number of nodes/polynomials
p = zeros(N,size(x,1));

r = x(:,1); s = x(:,2); 

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = ones(size(x,1),1); %1;
         q_i = ones(size(x,1),1); % 1;
     else
         p_i = jacobiPol(r,0,0,i);
         q_i = q_i.*(1-s)./2;
     end
     %Value for j
     j = nDeg-i;
     if j==0
        p_j = 1;
     else
        p_j = jacobiPol(s,2*i+1,0,j);
     end
     ncount = ncount+1;
     factor = sqrt( (2*i+1)*(i+j+1)/2 );
     p(ncount,:)    = ( p_i.*q_i.*p_j )*factor;
  end
end