function p = jacobiP(x,a,b,N)
% Jacoby polynomials
% p = jacobiP(x,a,b,n)

if N==0
    p = 1; return
end

if N==1
    p = (1/2)*( a-b + (2+a+b)*x ); return
end


apb = a+b;
apb_amb = (a+b)*(a-b);
pm2 = 1;
pm1 = (1/2)*( a-b + (2+a+b)*x );
for n = 2:N
  A = 2*n + apb; 
  B = n*(n+apb)*(A-2);
  p = ( (A-1)*( apb_amb + x*(A-2)*A)/(2*B) ).*pm1 ...
    - ( (n+a-1)*(n+b-1)*A/B ).*pm2;
  pm2 = pm1;
  pm1 = p;
end
    
