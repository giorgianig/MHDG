clear
clc

syms bn U1 U2 U3 U4
assume(bn,'real')
assume(U1,'real')
assume(U2,'real')
assume(U3,'real')
assume(U4,'real')

A= bn*([0,                                                                     1,                                       0,                                              0;...
            -(2/3)*U2^2/U1^2,                                   (4/3)*U2/U1,                             2/3,                                          2/3;...
    (2/3)*U2^3/U1^3-(5/3)*U2*U3/U1^2,    (5/3)*U3/U1-U2^2/U1^2,           (5/3)*U2/U1,                                       0;...
    -(5/3)*U2*U4/U1^2,                                         (5/3)*U4/U1,                               0,                             (5/3)*U2/U1]);



 [R,D] = eig(A);
 [Rt,Dt] = eig(A');
 
 disp([' transpose(Rt)*D*R  '])
 simplify(transpose(Rt)*D*R)
 
 disp([' transpose(Rt)*abs(D)*R  '])
  simplify(transpose(Rt)*abs(D)*R)
  
  disp([' R*abs(D)*inv(R) '])
 res = simplify(R*abs(D)*inv(R))
 
 disp((['fortran( R*abs(D)*inv(R) )']))
 fortran(res)
 
 disp([' A-R*D*inv(R)  '])
 simplify(A-R*D*inv(R))
 