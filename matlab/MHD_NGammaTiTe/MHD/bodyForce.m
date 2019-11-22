function [f, ftemp] = bodyForce(X)

global testcase Mesh diff_n diff_u diff_ei diff_ee diff_pari diff_pare tie axisym neq  epn pcourr
global Mref
D = diff_n;
mu = diff_u;
csii = diff_ei;
csie = diff_ee;
kpari = diff_pari;
kpare = diff_pare;



% allocation
f = zeros(size(X,1),neq);
ftemp = f;
% assign
x = X(:,1);
y = X(:,2);

xmax = Mesh.maxx;
xmin = Mesh.minx;
ymax = Mesh.maxy;
ymin = Mesh.miny;

xc = testcase.xc;
yc = testcase.yc;

a = testcase.wavex;
b = testcase.wavey;
a = a*2*pi;
b = b*2*pi;
xm = 0.5*(xmax+xmin);
ym = 0.5*(ymax+ymin);
xm = xc;
ym = yc;
%% bx = 0.1*Bx/|B| -- by = 0.1*By/|B|
switch testcase.n
    case 1
        if axisym,error('This is NOT an axisymmetric case!'),end
        for i=1:size(x,1)
            xx = x(i);
            yy = y(i);
            cx2 = cos(a*xx)^2;
            sx2 = sin(a*xx)^2;
            cy2 = cos(a*yy)^2;
            sy2 = sin(a*yy)^2;            
            t1 = cos(a*xx)*sin(a*yy);
            t2 = cos(a*yy)*sin(a*xx);
            t3 = ((xm - xx)^2 + (ym - yy)^2 + 1)^(1/2);
            t4 = ((xm - xx)^2 + (ym - yy)^2 + 1)^(3/2);            
            t5 = cx2*cy2;
            t6 = cos(a*xx)*cos(a*yy);
            t7 = sin(a*xx)*sin(a*yy);
            
        f(i,1) =  (2*D*a^2*t7 - 2*a*xm*t1*t3 + 2*a*xx*t1*t3 + 2*a*ym*t2*t3 - 2*a*yy*t2*t3 + D*a^2*xm^2*t7 + D*a^2*xx^2*t7 + D*a^2*ym^2*t7 + D*a^2*yy^2*t7 + D*a*xm*t1 - D*a*xx*t1 + D*a*ym*t2 - D*a*yy*t2 + a*xm*cos(a*xx)*cy2*sin(a*xx)*t3 - a*xx*cos(a*xx)*cy2*sin(a*xx)*t3 - a*ym*cx2*cos(a*yy)*sin(a*yy)*t3 + a*yy*cx2*cos(a*yy)*sin(a*yy)*t3 - a*xm*cos(a*xx)*sin(a*xx)*sy2*t3 + a*xx*cos(a*xx)*sin(a*xx)*sy2*t3 + a*ym*t2^2*sin(a*yy)*t3 - a*yy*t2^2*sin(a*yy)*t3 - 2*D*a^2*xm*ym*t6 + 2*D*a^2*xx*ym*t6 + 2*D*a^2*xm*yy*t6 - 2*D*a^2*xx*yy*t6 - 2*D*a^2*xm*xx*t7 - 2*D*a^2*ym*yy*t7)/((xm - xx)^2 + (ym - yy)^2 + 1);
        f(i,2) =  (a*(60*xm*t2*t3 - 60*xx*t2*t3 - 60*ym*t1*t3 + 60*yy*t1*t3 + 4*xm*t7*t3 - 4*xx*t7*t3 + 4*ym*t7*t3 - 4*yy*t7*t3 + 12*a*mu*t6 - 6*mu*xm*t2 + 6*mu*xx*t2 - 6*mu*ym*t1 + 6*mu*yy*t1 - 2*xm*cy2*sx2*t3 + 2*xx*cy2*sx2*t3 - 2*ym*cx2*sy2*t3 + 2*yy*cx2*sy2*t3 + 2*xm*sx2*sy2*t3 - 2*xx*sx2*sy2*t3 + 2*ym*sx2*sy2*t3 - 2*yy*sx2*sy2*t3 + 4*xm*t6*t3 - 4*xx*t6*t3 + 4*ym*t6*t3 - 4*yy*t6*t3 + 6*a*mu*xm^2*t6 + 6*a*mu*xx^2*t6 + 6*a*mu*ym^2*t6 + 6*a*mu*yy^2*t6 + 2*xm*cx2*cos(a*yy)^3*sin(a*xx)*t3 - 2*xx*cx2*cos(a*yy)^3*sin(a*xx)*t3 - 2*ym*cos(a*xx)^3*cy2*sin(a*yy)*t3 + 2*yy*cos(a*xx)^3*cy2*sin(a*yy)*t3 + 3*mu*ym*cos(a*xx)*cy2*sin(a*xx) + 3*mu*xm*cx2*cos(a*yy)*sin(a*yy) - 3*mu*xx*cx2*cos(a*yy)*sin(a*yy) - 3*mu*yy*cos(a*xx)*cy2*sin(a*xx) - 3*mu*ym*cos(a*xx)*sin(a*xx)*sy2 - 3*mu*xm*t2^2*sin(a*yy) + 3*mu*xx*t2^2*sin(a*yy) + 3*mu*yy*cos(a*xx)*sin(a*xx)*sy2 + 8*ym*cos(a*xx)*cy2*sin(a*xx)*t3 - 8*xm*cx2*cos(a*yy)*sin(a*yy)*t3 + 8*xx*cx2*cos(a*yy)*sin(a*yy)*t3 - 8*yy*cos(a*xx)*cy2*sin(a*xx)*t3 - 4*xm*cx2*t2*sy2*t3 + 4*xx*cx2*t2*sy2*t3 + 4*ym*cos(a*xx)*cy2*sx2*sin(a*yy)*t3 - 4*yy*cos(a*xx)*cy2*sx2*sin(a*yy)*t3 - 6*a*mu*xm*ym*t5 + 6*a*mu*xx*ym*t5 + 6*a*mu*xm*yy*t5 - 6*a*mu*xx*yy*t5 + 6*a*mu*xm*ym*cx2*sy2 + 6*a*mu*xm*ym*cy2*sx2 - 6*a*mu*xx*ym*cx2*sy2 - 6*a*mu*xx*ym*cy2*sx2 - 6*a*mu*xm*yy*cx2*sy2 - 6*a*mu*xm*yy*cy2*sx2 + 6*a*mu*xx*yy*cx2*sy2 + 6*a*mu*xx*yy*cy2*sx2 - 6*a*mu*xm*ym*sx2*sy2 + 6*a*mu*xx*ym*sx2*sy2 + 6*a*mu*xm*yy*sx2*sy2 - 6*a*mu*xx*yy*sx2*sy2 + 4*xm*cos(a*xx)*t2*sin(a*yy)*t3 - 4*xx*cos(a*xx)*t2*sin(a*yy)*t3 + 4*ym*cos(a*xx)*t2*sin(a*yy)*t3 - 4*yy*cos(a*xx)*t2*sin(a*yy)*t3 - 12*a*mu*xm*xx*t6 - 12*a*mu*ym*yy*t6 - 12*a*mu*xm*ym*t7 + 12*a*mu*xx*ym*t7 + 12*a*mu*xm*yy*t7 - 12*a*mu*xx*yy*t7 + 24*a*mu*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*xm^2*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*xx^2*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*ym^2*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*yy^2*cos(a*xx)*t2*sin(a*yy) - 24*a*mu*xm*xx*cos(a*xx)*t2*sin(a*yy) - 24*a*mu*ym*yy*cos(a*xx)*t2*sin(a*yy)))/(3*(xm^2 - 2*ym*yy - 2*xm*xx + xx^2 + ym^2 + yy^2 + 1));
        f(i,3) =   (csii*(a*xx - a*xm + 3*a^2*sin(2*a*xx) + 2*a*xm*cx2 - 2*a*xx*cx2 + a*xm*cy2 - a*xx*cy2 + 4*a^2*t1 + 40*a^2*t7 + 2*a^2*xm^2*sin(2*a*xx) + 2*a^2*xx^2*sin(2*a*xx) + a^2*ym^2*sin(2*a*xx) + a^2*yy^2*sin(2*a*xx) - 2*a*xm*t5 + 2*a*xx*t5 + 2*a^2*xm^2*t1 + 2*a^2*xx^2*t1 + 2*a^2*ym^2*t1 + 2*a^2*yy^2*t1 + 20*a^2*xm^2*t7 + 20*a^2*xx^2*t7 + 20*a^2*ym^2*t7 + 20*a^2*yy^2*t7 - 8*a^2*cos(a*xx)*cy2*sin(a*xx) + 2*a*ym*t6 - 2*a*yy*t6 + 20*a*xm*t1 - 20*a*xx*t1 + 20*a*ym*t2 - 20*a*yy*t2 - 2*a*xm*t7 + 2*a*xx*t7 - 4*a^2*xm*xx*sin(2*a*xx) + 2*a^2*xm*ym*sin(2*a*yy) - 2*a^2*xx*ym*sin(2*a*yy) - 2*a^2*ym*yy*sin(2*a*xx) - 2*a^2*xm*yy*sin(2*a*yy) + 2*a^2*xx*yy*sin(2*a*yy) - 4*a^2*xm^2*cos(a*xx)*cy2*sin(a*xx) - 4*a^2*xx^2*cos(a*xx)*cy2*sin(a*xx) - 4*a^2*ym^2*cos(a*xx)*cy2*sin(a*xx) - 4*a^2*yy^2*cos(a*xx)*cy2*sin(a*xx) - 40*a^2*xm*ym*t6 + 40*a^2*xx*ym*t6 + 40*a^2*xm*yy*t6 - 40*a^2*xx*yy*t6 - 4*a^2*xm*xx*t1 + 4*a^2*xm*ym*t2 - 4*a^2*xx*ym*t2 - 4*a^2*xm*yy*t2 + 4*a^2*xx*yy*t2 - 4*a^2*ym*yy*t1 - 40*a^2*xm*xx*t7 - 40*a^2*ym*yy*t7 + 8*a^2*xm*xx*cos(a*xx)*cy2*sin(a*xx) - 8*a^2*xm*ym*cx2*cos(a*yy)*sin(a*yy) + 8*a^2*xx*ym*cx2*cos(a*yy)*sin(a*yy) + 8*a^2*ym*yy*cos(a*xx)*cy2*sin(a*xx) + 8*a^2*xm*yy*cx2*cos(a*yy)*sin(a*yy) - 8*a^2*xx*yy*cx2*cos(a*yy)*sin(a*yy) + 2*a*ym*cos(a*xx)*t2*sin(a*yy) - 2*a*yy*cos(a*xx)*t2*sin(a*yy)))/(xm^2 - 2*ym*yy - 2*xm*xx + xx^2 + ym^2 + yy^2 + 1) - (kpari*(2*3^(1 - epn)*Mref*a^2*xx^3*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a^2*xx^3*t1*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*a^2*epn*xx^3*t5*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 2*3^(1 - epn)*Mref*a^2*xx^3*t5*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*xx^2*t7*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a*ym^2*t7*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a*yy^2*t7*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*a^2*epn*xx*ym^2*sx2*sy2*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 4*3^(1 - epn)*a^2*epn*xx*yy^2*sx2*sy2*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 2*3^(1 - epn)*Mref*a*xm*xx*t7*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*Mref*a*ym*yy*t7*((2*t1 - t5 + 40)/Mref)^epn + 8*3^(1 - epn)*a^2*epn*xx^3*cos(a*xx)^3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 4*3^(1 - epn)*Mref*a^2*xm*xx^2*t5*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xm^2*xx*t5*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xx*ym^2*t5*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xx*yy^2*t5*((2*t1 - t5 + 40)/Mref)^epn +...
            4*3^(1 - epn)*Mref*a^2*xm*xx^2*t1*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a^2*xm^2*xx*t1*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a^2*xx*ym^2*t1*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a^2*xx^2*ym*t2*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a^2*xx*yy^2*t1*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*Mref*a^2*xx^2*yy*t2*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a^2*xm*xx^2*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xm^2*xx*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xx*ym^2*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a^2*xx*yy^2*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a*xx^2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*ym^2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*yy^2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^epn - 8*3^(1 - epn)*a^2*epn*xm*xx^2*t5*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 4*3^(1 - epn)*a^2*epn*xm^2*xx*t5*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 2*3^(1 - epn)*Mref*a*xm*ym*t6*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*Mref*a*xx*ym*t6*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*xm*yy*t6*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a*xx*yy*t6*((2*t1 - t5 + 40)/Mref)^epn - 8*3^(1 - epn)*a^2*epn*xx*ym*yy*sx2*sy2*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 4*3^(1 - epn)*Mref*a^2*xx*ym*yy*t5*((2*t1 - t5 + 40)/Mref)^epn + 8*3^(1 - epn)*a^2*epn*xx^2*ym*cx2*cos(a*yy)^3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 16*3^(1 - epn)*a^2*epn*xm*xx^2*cos(a*xx)^3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*a^2*epn*xm^2*xx*cos(a*xx)^3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 8*3^(1 - epn)*a^2*epn*xx^2*yy*cx2*cos(a*yy)^3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 4*3^(1 - epn)*Mref*a^2*xm*xx*ym*t2*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a^2*xm*xx*yy*t2*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*Mref*a^2*xx*ym*yy*t1*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a^2*xx*ym*yy*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*xm*xx*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^epn - 2*3^(1 - epn)*Mref*a*xm*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn + 4*3^(1 - epn)*Mref*a*xx*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a*ym*yy*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^epn + 2*3^(1 - epn)*Mref*a*xm*yy*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn - 4*3^(1 - epn)*Mref*a*xx*yy*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn - 8*3^(1 - epn)*Mref*a^2*xx^2*ym*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn + 8*3^(1 - epn)*Mref*a^2*xx^2*yy*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn + 8*3^(1 - epn)*a^2*epn*xx^2*ym*cos(a*xx)^3*cos(a*yy)^3*t7*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 8*3^(1 - epn)*a^2*epn*xx^2*yy*cos(a*xx)^3*cos(a*yy)^3*t7*((2*t1 - t5 + 40)/Mref)^(epn - 1) - (4*3^(1 - epn)*Mref*a^2*epn*xx^3*cos(a*xx)^4*cy2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) - 8*3^(1 - epn)*a^2*epn*xx^2*ym*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*a^2*epn*xx^2*yy*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 8*3^(1 - epn)*a^2*epn*xm*xx*ym*cx2*cos(a*yy)^3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*a^2*epn*xm*xx*yy*cx2*cos(a*yy)^3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*Mref*a^2*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn - 8*3^(1 - epn)*Mref*a^2*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^epn - 8*3^(1 - epn)*a^2*epn*xm*xx*ym*cos(a*xx)^3*cos(a*yy)^3*t7*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*a^2*epn*xm*xx*yy*cos(a*xx)^3*cos(a*yy)^3*t7*((2*t1 - t5 + 40)/Mref)^(epn - 1) + 8*3^(1 - epn)*a^2*epn*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - 8*3^(1 - epn)*a^2*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)^(epn - 1) - (4*3^(1 - epn)*Mref*a^2*epn*xx*ym^2*cx2*cos(a*yy)^4*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xm*xx^2*cos(a*xx)^4*cy2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) - (4*3^(1 - epn)*Mref*a^2*epn*xm^2*xx*cos(a*xx)^4*cy2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) - (4*3^(1 - epn)*Mref*a^2*epn*xx*yy^2*cx2*cos(a*yy)^4*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xx*ym*yy*cx2*cos(a*yy)^4*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xx*ym^2*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xx^2*ym*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xx*yy^2*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) - (8*3^(1 - epn)*Mref*a^2*epn*xx^2*yy*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) - (8*3^(1 - epn)*Mref*a^2*epn*xm*xx*ym*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - t5 + 40) + (8*3^(1 - epn)*Mref*a^2*epn*xm*xx*yy*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/(2*t1 - ...
            t5 + 40) - (16*3^(1 - epn)*Mref*a^2*epn*xx*ym*yy*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)^epn)/...
            (2*t1 - t5 + 40)))/(9*Mref^2*xx*((xm - xx)^2 + (ym - yy)^2 + 1)) + Mref*pcourr*t6*((2*a*sin(a*xx)*(xm - xx)*...
            (10*cos(a*yy) + sin(a*xx) + 2*sin(a*yy) - 2*cy2*sin(a*xx)))/(3*Mref*t3) + (4*a*cos(a*xx)*(ym - yy)*(cos(a*yy) - 5*sin(a*yy) + ...
            t2*sin(a*yy)))/(3*Mref*t3)) - (a*t6*(ym - yy)*(5*cx2*sy2 - 5*sx2*sy2 +...
            100*t1 - 10*t7 - cos(a*xx)^3*cy2*sin(a*yy) + 4*cos(a*xx)*cy2*sin(a*xx) + 2*cos(a*xx)*cy2*sx2*sin(a*yy)))/(3*((xm - xx)^2 + ...
            (ym - yy)^2 + 1)^(1/2)) + (a*cos(a*xx)*cy2*(xm - xx)*(10*cos(a*xx) + 100*sin(a*xx) + 2*cx2*sin(a*xx) + 4*cx2*sin(a*yy) - ...
            3*t5*sin(a*xx) + 10*cos(a*xx)*t7))/(3*t3) + (3^(1/2)*(t7 + 2)^2*(2*t1 -...
            t5 + 2*t2 + 20))/(2*tie*(t2 - 10)*(-(2*t2 - 20)/Mref)^(1/2)) -...
            (a*t1*(t7 + 2)*(xm - xx)*(5*t1 - t5 + 100))/(3*((xm - xx)^2 + ...
            (ym - yy)^2 + 1)^(1/2)) + (a*t2*(t7 + 2)*(ym - yy)*(5*t1 - t5 + 100))/(3*t3) - ...
            (t6*(2*xm - 2*xx)*(t7 + 2)*(ym - yy)*(5*t1 - t5 + 100))/(6*t4) + ...
            (t6*(2*ym - 2*yy)*(t7 + 2)*(xm - xx)*(5*t1 - t5 + 100))/(6*t4);
        f(i,4) = (csie*(a*ym - a*yy - 3*a^2*sin(2*a*yy) - a*ym*cx2 + a*yy*cx2 - 2*a*ym*cy2 + 2*a*yy*cy2 - 4*a^2*t2 + 20*a^2*t7 - a^2*xm^2*sin(2*a*yy) - a^2*xx^2*sin(2*a*yy) - 2*a^2*ym^2*sin(2*a*yy) - 2*a^2*yy^2*sin(2*a*yy) + 2*a*ym*t5 - 2*a*yy*t5 - 2*a^2*xm^2*t2 - 2*a^2*xx^2*t2 - 2*a^2*ym^2*t2 - 2*a^2*yy^2*t2 + 10*a^2*xm^2*t7 + 10*a^2*xx^2*t7 + 10*a^2*ym^2*t7 + 10*a^2*yy^2*t7 + 8*a^2*cx2*cos(a*yy)*sin(a*yy) - 2*a*xm*t6 + 2*a*xx*t6 + 10*a*xm*t1 - 10*a*xx*t1 + 10*a*ym*t2 - 10*a*yy*t2 + 2*a*ym*t7 - 2*a*yy*t7 - 2*a^2*xm*ym*sin(2*a*xx) + 2*a^2*xx*ym*sin(2*a*xx) + 2*a^2*xm*xx*sin(2*a*yy) + 2*a^2*xm*yy*sin(2*a*xx) - 2*a^2*xx*yy*sin(2*a*xx) + 4*a^2*ym*yy*sin(2*a*yy) + 4*a^2*xm^2*cx2*cos(a*yy)*sin(a*yy) + 4*a^2*xx^2*cx2*cos(a*yy)*sin(a*yy) + 4*a^2*ym^2*cx2*cos(a*yy)*sin(a*yy) + 4*a^2*yy^2*cx2*cos(a*yy)*sin(a*yy) - 20*a^2*xm*ym*t6 + 20*a^2*xx*ym*t6 + 20*a^2*xm*yy*t6 - 20*a^2*xx*yy*t6 + 4*a^2*xm*xx*t2 - 4*a^2*xm*ym*t1 + 4*a^2*xx*ym*t1 + 4*a^2*xm*yy*t1 - 4*a^2*xx*yy*t1 + 4*a^2*ym*yy*t2 - 20*a^2*xm*xx*t7 - 20*a^2*ym*yy*t7 + 8*a^2*xm*ym*cos(a*xx)*cy2*sin(a*xx) - 8*a^2*xx*ym*cos(a*xx)*cy2*sin(a*xx) - 8*a^2*xm*xx*cx2*cos(a*yy)*sin(a*yy) - 8*a^2*xm*yy*cos(a*xx)*cy2*sin(a*xx) + 8*a^2*xx*yy*cos(a*xx)*cy2*sin(a*xx) - 8*a^2*ym*yy*cx2*cos(a*yy)*sin(a*yy) - 2*a*xm*cos(a*xx)*t2*sin(a*yy) + 2*a*xx*cos(a*xx)*t2*sin(a*yy)))/(xm^2 - 2*ym*yy - 2*xm*xx + xx^2 + ym^2 + yy^2 + 1) - Mref*pcourr*t6*((2*a*sin(a*xx)*(xm - xx)*(10*cos(a*yy) + sin(a*xx) + 2*sin(a*yy) - 2*cy2*sin(a*xx)))/(3*Mref*t3) + (4*a*cos(a*xx)*(ym - yy)*(cos(a*yy) - 5*sin(a*yy) + t2*sin(a*yy)))/(3*Mref*t3)) + (10*a*cx2*cos(a*yy)*(ym - yy)*(cos(a*yy) - 5*sin(a*yy) + t2*sin(a*yy)))/(3*t3) - (3^(1/2)*(t7 + 2)^2*(2*t1 - t5 + 2*t2 + 20))/(2*tie*(t2 - 10)*(-(2*t2 - 20)/Mref)^(1/2)) + (5*a*cos(a*xx)*t2*(xm - xx)*(10*cos(a*yy) + sin(a*xx) + 2*sin(a*yy) - 2*cy2*sin(a*xx)))/(3*t3) + (5*t6*(t2 - 10)*(2*xm - 2*xx)*(t7 + 2)*(ym - yy))/(6*t4) - (5*t6*(t2 - 10)*(2*ym - 2*yy)*(t7 + 2)*(xm - xx))/(6*t4) + (5*a*t1*(t2 - 10)*(t7 + 2)*(xm - xx))/(3*t3) - (5*a*t2*(t2 - 10)*(t7 + 2)*(ym - yy))/(3*t3) + (2/3^(epn + 1)*a*kpare*(-(2*t2 - 20)/Mref)^epn*(10*xx^2*t6 - 10*ym^2*t6 - 10*yy^2*t6 - xx^2*cos(a*xx)*cy2*sin(a*xx) + ym^2*cos(a*xx)*cy2*sin(a*xx) + yy^2*cos(a*xx)*cy2*sin(a*xx) - 10*xm*xx*t6 + 20*ym*yy*t6 - a*xx^3*cy2*sx2 - 10*xm*ym*t7 + 20*xx*ym*t7 + 10*xm*yy*t7 - 20*xx*yy*t7 + 10*a*xx^3*t2 + a*epn*xx^3*sx2*sy2 + 2*a*xm*xx^2*cy2*sx2 - a*xm^2*xx*cy2*sx2 - a*xx*ym^2*cy2*sx2 - a*xx*yy^2*cy2*sx2 - 20*a*xm*xx^2*t2 + 10*a*xm^2*xx*t2 + 10*a*xx*ym^2*t2 + 20*a*xx^2*ym*t1 + 10*a*xx*yy^2*t2 - 20*a*xx^2*yy*t1 + xm*xx*cos(a*xx)*cy2*sin(a*xx) - 2*ym*yy*cos(a*xx)*cy2*sin(a*xx) + xm*ym*t2^2*sin(a*yy) - 2*xx*ym*t2^2*sin(a*yy) - xm*yy*t2^2*sin(a*yy) + 2*xx*yy*t2^2*sin(a*yy) + 2*a*xx*ym*yy*cy2*sx2 + a*epn*xx*ym^2*t5 + a*epn*xx*yy^2*t5 - 20*a*xm*xx*ym*t1 + 20*a*xm*xx*yy*t1 - 20*a*xx*ym*yy*t2 - 2*a*epn*xm*xx^2*sx2*sy2 + a*epn*xm^2*xx*sx2*sy2 - 2*a*epn*xx*ym*yy*t5 - 2*a*xx^2*ym*cos(a*xx)*t2*sin(a*yy) + 2*a*xx^2*yy*cos(a*xx)*t2*sin(a*yy) - 2*a*epn*xx^2*ym*cos(a*xx)*t2*sin(a*yy) + 2*a*epn*xx^2*yy*cos(a*xx)*t2*sin(a*yy) + 2*a*xm*xx*ym*cos(a*xx)*t2*sin(a*yy) - 2*a*xm*xx*yy*cos(a*xx)*t2*sin(a*yy) + 2*a*epn*xm*xx*ym*cos(a*xx)*t2*sin(a*yy) - 2*a*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)))/(Mref*xx*(t2 - 10)*(xm^2 - 2*ym*yy - 2*xm*xx + xx^2 + ym^2 + yy^2 + 1));


        end
    case 2
        %% bx = 1/30; by = 1/30; n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
        f = bf2(X);
        
    case 3
        %%  bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n = 1+x^2+y^2, u = x^2+y^2; E = 10+x^2+y^2;
        
        f = bf3(X);
    case 4
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        f = bf4(X);
    case 5 % test for Bohm bc
        %% bx = 1/30 ; by = 0; n=1.1-sin(xr);u=sin(xr); Ei = 2-sin(xr); Ee = 1.4-sin(xr);
        f = bf5(X);
    case 6 % test for Bohm bc
        f = bf6(X);       
    case 25
        % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        f = bf25(X);
    case 26
        % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20-sin(a*xr)*cos(b*yr);
        f = bf26(X);
        
    case 30
        %% Circular field centered in [0, 0], n: density spot, u=0
    case 31
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
    case 32
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
        
        %% West
    case 50
        % WEST
        f(:,1) = 0;
        f(:,2) = 0;
    case 51
        % WEST
        error('You should not be here!')
    case 52
        % WEST
        error('You should not be here!')
    case 53
        % WEST
        error('You should not be here!')
    case 54
        % WEST
        error('You should not be here!')
    case 55
        % WEST
        error('You should not be here!')
        %% Circular case with limiter
    case 60
        %
        f(:,1) = 0;
        f(:,2) = 0;
        %% Circular case with limiter all Dirichlet bc
    case 69
        %
        f(:,1) = 0;
        f(:,2) = 0;        
    case 100
        if axisym,error('This is NOT an axisymmetric case!'),end
        %% Bx = 1, By = 0, n = 1, u = 0, E= 1;
    otherwise
        error('case not valid')
        
end
