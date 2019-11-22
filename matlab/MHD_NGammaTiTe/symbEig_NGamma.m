clear
clc

syms bn U1 U2 U3 U4 k
assume(bn,'real')
assume(k,'real')
assume(U1,'real')
assume(U2,'real')


A= bn*([0,                                                                     1;...
            k-U2^2/U1^2,                                   2*U2/U1]);



 [R,D] = eig(A);
%  [Rt,Dt] = eig(A');
 
%  disp([' transpose(Rt)*D*R  '])
%  simplify(transpose(Rt)*D*R)
%  
%  disp([' transpose(Rt)*abs(D)*R  '])
%   simplify(transpose(Rt)*abs(D)*R)
%   
%   disp([' R*abs(D)*inv(R) '])
%  simplify(R*abs(D)*inv(R))
%  
%  disp([' A-R*D*inv(R)  '])
%  simplify(A-R*D*inv(R))