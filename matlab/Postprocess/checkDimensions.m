clear

u = symunit;

L0 =  1.901e-3*u.m;
t0 = 7.28e6^-1*u.s;
n0 = 1e19*u.m^-3;
u0 = L0/t0;
m0 = 3.35e-27*u.kg;
e  = 1.60217662e-19*u.A*u.s;
Jp0  = u.A*u.m^-2;
phi0 = u.kg*u.m^2/u.A/u.s^3;
kb = 1.38064852e-23*u.m^2*u.kg*u.s^-2*u.K^-1;
B0 = u.kg*u.A^-1*u.s^-2;
W0 = u.A*u.s*u.kg^-1;
etap = u.m^3*u.kg*u.s^-3*u.A^-2;
T0 = 5.802261027303315e+05*u.K;
Tev = T0*kb/e;


T0*W0/u0^2

n0*t0*m0*W0^2

T0*t0/(L0^2*n0*m0*W0)

1/W0/t0

n0*u0^3*t0*L0^2/W0/T0


 unitInfo((2*Tev*e/m0)/(L0/t0)^2)
 
 unitInfo(phi0/B0*t0/L0^2)
 
 unitInfo( Tev*t0/(L0^2*B0))

 
