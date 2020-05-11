clear
home

test =2;
syms xx yy xm ym tt br bz bt Br Bz Bt Bp BB divb divb_cyl b
syms n u a theta
syms PaGradN PeGradNx PeGradNy PeGradNt
syms DivPeGradN D k
syms PaGradNu PeGradNux PeGradNuy PeGradNut coeffconv
syms DivPeGradNu mu
syms f1 f2
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
    case(2)
        Br =  (yy-ym)/xx;
        Bz = -(xx-xm)/xx;
        Bt = 1;
        axisym = 1;
        n = 2+sin(a*xx)*sin(a*yy);
        u = cos(a*xx)*cos(a*yy);
    case(3)
        Br =  0*xx;
        Bz =  1;
        Bt =  1;
        axisym = 0;
        n = 1+xx;
        u = yy;
    case(4)
        Br =  0*xx;
        Bz =  1;
        Bt =  1;
        axisym = 0;
        n = 1;
        u = xx;
    case(6)
        Br =  0*xx;
        Bz =  1;
        Bt =  1;
        axisym = 0;
        n = 1+xx;
        u = xx^3;
    case(7)
        Br =  0*xx;
        Bz =  1;
        Bt =  1;
        axisym = 0;
        n = 1+xx+yy;
        u = xx^2+yy^2;
        
    case(8)
        Br =  (yy-ym);
        Bz = -(xx-xm);
        Bt = 1;
        axisym = 0;
        n = 1+xx;
        u = 1+yy;
        
    case(9)
        Br =  0*xx;
        Bz =1;
        Bt = 1;
        axisym = 0;
        n = 1+xx;
        u = 1+yy;
    case(10)
        Br =  (yy-ym);
        Bz = -(xx-xm);
        Bt = 1;
        axisym = 0;
        n = 1;
        u = 1;
    case 100
        Bt = 1;
        Br = 2*pi/theta*(yy-ym);
        Bz = 2*pi/theta*(-xx+xm);
        n = 2+sin(a* xx)*sin(a* yy)*sin(b*tt);
        u = cos(a* xx)*cos(a* yy)*cos(b*tt);
        axisym = 0;
    case 1000
        Bt = 1;
        Br =  (yy-ym)/xx;
        Bz = -(xx-xm)/xx;
        n = 1;
        u = 1;
        axisym = 1;
    case 1001
        Bt = 1;
        Br = 1/xx;
        Bz = 0;
        n = 1;
        u = 1;
        axisym = 1;
    case 1002
        Bt = 1;
        Br = 0*xx;
        Bz = 0*xx;
        n = 1;
        u = 1;
        axisym = 1;       
    case 1003
        Bt = 0;
        Br = 1/xx;
        Bz = 0;
        n = 1;
        u = 1;
        axisym = 1;         
    case(2000)
        Br =  (yy-ym)/xx;
        Bz = -(xx-xm)/xx;
        Bt = 1;
        axisym = 1;
        n = 1;
        u = cos(a*xx)*cos(a*yy);        
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
    f1 =  (1/xx*diff(xx*n*u*br, xx)+diff(n*u*bz, yy)+1/xx*diff(n*u*bt, tt))-DivPeGradN;
else
    f1 =  diff(n*u*br, xx)+diff(n*u*bz, yy)+diff(n*u*bt, tt)-DivPeGradN;
end
%% second equation
% Parallel gradient of nn (scalar)
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
if axisym
    f2 =  (1/xx*diff(xx*n*u^2*br, xx)+diff(n*u^2*bz, yy)+1/xx*diff(n*u^2*bt, tt)+ k*br*(diff(n, xx))+k*bz*(diff(n, yy))+k*bt*(diff(n, tt))/xx) -DivPeGradNu;
else
    f2 =  diff(n*u^2*br, xx)+diff(n*u^2*bz, yy)+diff(n*u^2*bt, tt)+k*br*(diff(n, xx))+k*bz*(diff(n, yy))+k*bt*(diff(n, tt)) -DivPeGradNu;
end
disp([' f1:  ',  char(simplify(f1))])
disp([' f2:  ',  char(simplify(f2))])