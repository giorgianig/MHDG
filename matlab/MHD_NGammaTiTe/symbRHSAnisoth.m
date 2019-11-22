clear
home

test =2;
syms xx yy xm ym tt br bz bt Br Bz Bt Bp BB divb divb_cyl
syms n u a Ei Ee pio pel Ti Te
syms PaGradN PeGradNx PeGradNy PeGradNt
syms DivPeGradN D k
syms PaGradNu PeGradNux PeGradNuy PeGradNut
syms DivPeGradNu mu csii csie kpari kpare pcourr Mref epn tie
syms PaGradNEi PeGradNEix PeGradNEiy PeGradNEit DivPeGradNEi DivGradParTi GradParPe
syms PaGradNEe PeGradNEex PeGradNEey PeGradNEet DivPeGradNEe DivGradParTe
syms f1 f2 f3 f4
assume(br,'real')
assume(bz,'real')
assume(bt,'real')
assume(Br,'real')
assume(Bz,'real')
assume(Bt,'real')
assume(Bp,'real')
assume(BB,'real')
assume(divb,'real')
assume(divb_cyl,'real')

switch test
    case(1)
        Br =  (yy-ym);
        Bz = -(xx-xm);
        Bt = 1;
        axisym = 0;
        n = 2+sin(a*xx)*sin(a*yy);
        u = cos(a*xx)*cos(a*yy);
        Ei = 20+cos(a*xx)*sin(a*yy);
        Ee = 10-sin(a*xx)*cos(a*yy);
        
    case(2)
        Br =  (yy-ym)/xx;
        Bz = -(xx-xm)/xx;
        Bt = 1;       
        axisym = 1;
        n = 2+sin(a*xx)*sin(a*yy);
        u = cos(a*xx)*cos(a*yy);      
        Ei = 20+cos(a*xx)*sin(a*yy);
        Ee = 10-sin(a*xx)*cos(a*yy);    
        
    otherwise
        error('wrong case')
end

% check div(B) = 0
if axisym
    divB = 1/xx*diff(xx*Br,xx) + diff(Bz,yy) + 1/xx*diff(Bt,tt);
else
    divB = diff(Br,xx) + diff(Bz,yy) + diff(Bt,tt);
end

disp([' CHECK--->div(B):  ', char(simplify(divB))])

Bp = sqrt(Br^2+Bz^2);
BB = sqrt(Bp^2+Bt^2);

br = Br/BB;
bz = Bz/BB;
bt = Bt/BB;

divb = diff(br,xx) + diff(bz,yy) + diff(bt,tt);
divb_cyl = 1/xx*diff(br*xx,xx)+diff(bz,yy)+1/xx*diff(bt,tt);

disp([' br:  ',  char(simplify(br))])
disp([' bz:  ', char(simplify(bz))])
disp([' bt:  ', char(simplify(bt))])
disp([' divb:  ', char(simplify(divb))])
disp([' divb cyl:  ',char(simplify(divb_cyl))])


%% Temperature and pressure
pio = 2*n*(Ei-(1/2)*u^2)/(3*Mref);
pel = 2*n*Ee/(3*Mref);
Ti = pio/n;
Te = pel/n;

%% first equation
% Parallel gradient of n (scalar)
if axisym
    PaGradN = br*(diff(n, xx))+bz*(diff(n, yy))+bt*(diff(n, tt))/xx;
else
    PaGradN = br*(diff(n, xx))+bz*(diff(n, yy))+bt*(diff(n, tt));
end
% Perpendicular gradient of n (vector)
PeGradNx = diff(n, xx)-br*PaGradN;
PeGradNy = diff(n, yy)-bz*PaGradN;
if axisym
    PeGradNt  = diff(n, tt)/xx-bt*PaGradN;
else
    PeGradNt  = diff(n, tt)-bt*PaGradN;
end
% Divergence of the perpendicular gradient of n multiplied by diffusion
if axisym
    DivPeGradN = D*(1/xx*diff(xx*PeGradNx, xx)+diff(PeGradNy, yy))+1/xx*diff(PeGradNt, tt);
else
    DivPeGradN = D*(diff(PeGradNx, xx)+diff(PeGradNy, yy))+diff(PeGradNt, tt);
end

if axisym
    f1 =  1/xx*diff(xx*n*u*br, xx)+diff(n*u*bz, yy)+1/xx*diff(n*u*bt, tt)-DivPeGradN;
else
    f1 =  diff(n*u*br, xx)+diff(n*u*bz, yy)+diff(n*u*bt, tt)-DivPeGradN;
end


%% second equation
% Parallel gradient of nu (scalar)
if axisym
    PaGradNu = br*(diff(n*u, xx))+bz*(diff(n*u, yy))+bt*(diff(n*u, tt))/xx;
else
    PaGradNu = br*(diff(n*u, xx))+bz*(diff(n*u, yy))+bt*(diff(n*u, tt));
end
% Perpendicular gradient of n (vector)
PeGradNux = diff(n*u, xx)-br*PaGradNu;
PeGradNuy = diff(n*u, yy)-bz*PaGradNu;
if axisym
    PeGradNut  = diff(n*u, tt)/xx-bt*PaGradNu;
else
    PeGradNut  = diff(n*u, tt)-bt*PaGradNu;
end

% Divergence of the perpendicular gradient of n multiplied by diffusion
if axisym
    DivPeGradNu = mu*(1/xx*diff(xx*PeGradNux, xx)+diff(PeGradNuy, yy))+1/xx*diff(PeGradNut, tt);
else
    DivPeGradNu = mu*(diff(PeGradNux, xx)+diff(PeGradNuy, yy))+diff(PeGradNut, tt);
end

% Parallel gradient of total pressure
if axisym
    GradParPress = br*(diff(pio+pel, xx))+bz*(diff(pio+pel, yy))+bt*(diff(pio+pel, tt))/xx;
else
    GradParPress = br*(diff(pio+pel, xx))+bz*(diff(pio+pel, yy))+bt*(diff(pio+pel, tt));
end
    

if axisym
    f2 =  1/xx*diff(xx*n*u^2*br, xx)+diff(n*u^2*bz, yy)+1/xx*diff(n*u^2*bt, tt)+ Mref*GradParPress -DivPeGradNu;
else
    f2 =  diff(n*u^2*br, xx)+diff(n*u^2*bz, yy)+diff(n*u^2*bt, tt)+Mref*GradParPress -DivPeGradNu;
end



%% Third equation
% Parallel gradient of nEi (scalar)
if axisym
    PaGradNEi = br*(diff(n*Ei, xx))+bz*(diff(n*Ei, yy))+bt*(diff(n*Ei, tt))/xx;
else
    PaGradNEi = br*(diff(n*Ei, xx))+bz*(diff(n*Ei, yy))+bt*(diff(n*Ei, tt));
end

% Perpendicular gradient of nEi (vector)
PeGradNEix = diff(n*Ei, xx)-br*PaGradNEi;
PeGradNEiy = diff(n*Ei, yy)-bz*PaGradNEi;
if axisym
    PeGradNEit  = diff(n*Ei, tt)/xx-bt*PaGradNEi;
else
    PeGradNEit  = diff(n*Ei, tt)-bt*PaGradNEi;
end

% Divergence of the perpendicular gradient of nEi multiplied by diffusion
if axisym
    DivPeGradNEi = csii*(1/xx*diff(xx*PeGradNEix, xx)+diff(PeGradNEiy, yy))+1/xx*diff(PeGradNEit, tt);
else
    DivPeGradNEi = csii*(diff(PeGradNEix, xx)+diff(PeGradNEiy, yy))+diff(PeGradNEit, tt);
end

% Parallel gradient of the ionic temperature
if axisym
    GradParTempi = br*(diff(Ti, xx))+bz*(diff(Ti, yy))+bt*(diff(Ti, tt))/xx;
else
    GradParTempi = br*(diff(Ti, xx))+bz*(diff(Ti, yy))+bt*(diff(Ti, tt));
end

% Divergence of k//T^(epn)parallel gradient of Ti b
if axisym
    DivGradParTi = kpari*(1/xx*diff(xx*Ti^epn*GradParTempi*br, xx)+diff(Ti^epn*GradParTempi*bz, yy)+1/xx*diff(Ti^epn*GradParTempi*bt, tt));
else
    DivGradParTi = kpari*(diff(Ti^epn*GradParTempi*br, xx)+diff(Ti^epn*GradParTempi*bz, yy)+diff(Ti^epn*GradParTempi*bt, tt));
end

%% Parallel gradient of electronic pressure
if axisym
    GradParPe = br*(diff(pel, xx))+bz*(diff(pel, yy))+bt*(diff(pel, tt))/xx;
else
    GradParPe = br*(diff(pel, xx))+bz*(diff(pel, yy))+bt*(diff(pel, tt));
end

if axisym
    f3 =  1/xx*diff(xx*(n*Ei+Mref*pio)*u*br, xx)+diff((n*Ei+Mref*pio)*u*bz, yy)+1/xx*diff((n*Ei+Mref*pio)*u*bt, tt)...
           -DivPeGradNEi-DivGradParTi+pcourr*Mref*u*GradParPe+n^2*(Te-Ti)/(tie*Te^(3/2));
else
    f3 =  diff((n*Ei+Mref*pio)*u*br, xx)+diff((n*Ei+Mref*pio)*u*bz, yy)+diff((n*Ei+Mref*pio)*u*bt, tt)...
           -DivPeGradNEi-DivGradParTi+pcourr*Mref*u*GradParPe+n^2*(Te-Ti)/(tie*Te^(3/2));
end





%% Fourth equation
% Parallel gradient of nEe (scalar)
if axisym
    PaGradNEe = br*(diff(n*Ee, xx))+bz*(diff(n*Ee, yy))+bt*(diff(n*Ee, tt))/xx;
else
    PaGradNEe = br*(diff(n*Ee, xx))+bz*(diff(n*Ee, yy))+bt*(diff(n*Ee, tt));
end

% Perpendicular gradient of nEe (vector)
PeGradNEex = diff(n*Ee, xx)-br*PaGradNEe;
PeGradNEey = diff(n*Ee, yy)-bz*PaGradNEe;
if axisym
    PeGradNEet  = diff(n*Ee, tt)/xx-bt*PaGradNEe;
else
    PeGradNEet  = diff(n*Ee, tt)-bt*PaGradNEe;
end

% Divergence of the perpendicular gradient of nEi multiplied by diffusion
if axisym
    DivPeGradNEe = csie*(1/xx*diff(xx*PeGradNEex, xx)+diff(PeGradNEey, yy))+1/xx*diff(PeGradNEet, tt);
else
    DivPeGradNEe = csie*(diff(PeGradNEex, xx)+diff(PeGradNEey, yy))+diff(PeGradNEet, tt);
end

% Parallel gradient of the electronic temperature
if axisym
    GradParTempe = br*(diff(Te, xx))+bz*(diff(Te, yy))+bt*(diff(Te, tt))/xx;
else
    GradParTempe = br*(diff(Te, xx))+bz*(diff(Te, yy))+bt*(diff(Te, tt));
end

% Divergence of k//T^(epn)parallel gradient of Te b
if axisym
    DivGradParTe = kpare*(1/xx*diff(xx*Te^epn*GradParTempe*br, xx)+diff(Te^epn*GradParTempe*bz, yy)+1/xx*diff(Te^epn*GradParTempe*bt, tt));
else
    DivGradParTe = kpare*(diff(Te^epn*GradParTempe*br, xx)+diff(Te^epn*GradParTempe*bz, yy)+diff(Te^epn*GradParTempe*bt, tt));
end

if axisym
    f4 =  1/xx*diff(xx*(n*Ee+Mref*pel)*u*br, xx)+diff((n*Ee+Mref*pel)*u*bz, yy)+1/xx*diff((n*Ee+Mref*pel)*u*bt, tt)...
           -DivPeGradNEe-DivGradParTe-pcourr*Mref*u*GradParPe-n^2*(Te-Ti)/(tie*Te^(3/2));
else
     f4 =  diff((n*Ee+Mref*pel)*u*br, xx)+diff((n*Ee+Mref*pel)*u*bz, yy)+diff((n*Ee+Mref*pel)*u*bt, tt)...
           -DivPeGradNEe-DivGradParTe-pcourr*Mref*u*GradParPe-n^2*(Te-Ti)/(tie*Te^(3/2));
end

disp([' f1:  ',  char(simplify(f1))])
disp([' f2:  ',  char(simplify(f2))])
disp([' f3:  ',  char(simplify(f3))])
disp([' f4:  ',  char(simplify(f4))])



