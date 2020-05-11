syms x y r B0 Br Bz B Bt R0 t q br bz bt br_s bz_s

assume(br,'real')
assume(bz,'real')
assume(bt,'real')
assume(Br,'real')
assume(Bz,'real')
assume(Bt,'real')
assume(B,'real')
assume(x,'positive')
assume(y,'real')
assume(t,'real')
assume(r,'positive')
assume(q,'positive')
assume(B0,'positive')
assume(R0,'positive')
r  = sqrt((x-R0)^2+y^2);
Br = -B0*y/(x*q*sqrt(1- (r/R0)^2 ) );
Bz = B0*(x-R0)/(x*q*sqrt(1- (r/R0)^2 ) );
Bt = B0*R0/x;

B = sqrt(Br^2+Bz^2+Bt^2);
br = simplify(Br/B);
bz = simplify(Bz/B);
bt = simplify(Bt/B);

simplify(1/x*diff(x*br,x)+diff(bz,y)+1/x*diff(bt,t),'Steps',100,'IgnoreanalyticConstraints',true)



br_s = -y/sqrt(R0^2*q^2+(1-q^2)*r^2);
bz_s = (x-R0)/sqrt(R0^2*q^2+(1-q^2)*r^2);


simplify(1/x*diff(x*br_s,x)+diff(bz_s,y),'Steps',100,'IgnoreanalyticConstraints',true)