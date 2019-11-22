clear 
home

syms xr yr tr R0 B0 r Br Bz Bt q
r  = sqrt((xr-R0)^2+yr^2);
Br = -B0*yr/(xr*q*sqrt(1- (r/R0)^2 ) );
Bz = B0*(xr-R0)/(xr*q*sqrt(1- (r/R0)^2 ) );
Bt = B0*R0/xr;




Bp = sqrt(Br^2+Bz^2);
BB = sqrt(Bp^2+Bt^2);
br = Br/BB;
bz = Bz/BB;
bt = Bt/BB;


divB = simplify(1/xr*diff(xr*br,xr) + diff(bz,yr) + 1/xr*diff(bt,tr))