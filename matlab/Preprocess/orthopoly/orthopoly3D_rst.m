function p = orthopoly3D_rst(x,n)
% p = orthopoly3D_rst(x,n)
% Computes the ortogonal base of 3D polynomials of degree less 
% or equal to n at the point x=(r,s,t) in [-1,1]^3
%

N = (n+1)*(n+2)*(n+3)/6 ;%number of nodes/polynomials
p = zeros(N,1);

r = x(1); s = x(2); t = x(3);

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = 1;  q_i = 1;  
     else
         p_i = jacobiPol(r,0,0,i);     q_i = q_i*(1-s)/2;
     end
     %Loop increasing j
     for j = 0:(nDeg-i)
        if j==0
           p_j = 1;  q_j = ((1-t)/2)^i;
        else
           p_j = jacobiPol(s,2*i+1,0,j);    q_j = q_j*(1-t)/2;
        end
        %Value for k
        k = nDeg-(i+j);
        if k==0
           p_k = 1;
        else
           p_k = jacobiPol(t,2*(i+j)+2,0,k);
        end
        ncount = ncount+1;
        factor = sqrt( (2*i+1)*(i+j+1)*(2*(i+j+k)+3)/4 );
        p(ncount) = ( p_i*q_i*p_j*q_j*p_k )*factor;
    end
  end
end