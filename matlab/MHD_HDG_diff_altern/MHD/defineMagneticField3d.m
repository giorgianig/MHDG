function [b,db,dvx,dvy] = defineMagneticField3d(X,t)

global testcase Magnetic Mesh axisym driftvel theta
% allocation
nd = 3;
N2d = size(X,1);
N1d = size(t,1);
b = zeros(N2d*N1d,nd);
xmax = max(Mesh.X(:,1));
xmin = min(Mesh.X(:,1));
ymax = max(Mesh.X(:,2));
ymin = min(Mesh.X(:,2));
xm = 0.5*(xmax+xmin);
ym = 0.5*(ymax+ymin);

x = X(:,1);
y = X(:,2);

br = zeros(N2d,N1d);
bz = br;
bt = br;

db = br; % dbx_dx, dby_dy

xc = testcase.xc;
yc = testcase.yc;

% pitch angle
p = 85/180*pi;
for i=1:N2d
    for j=1:N1d
        xx = x(i);
        yy = y(i);
        tt = t(j);
        switch testcase.n
            case 1
                Br =  (yy-yc);
                Bz = -(xx-xc);
                Bt = 1;
                Bp = sqrt(Br^2+Bz^2);
                divb = 0;
            case 2
                Bt = 1;
                Br =  (yy-ym)/xx;
                Bz = -(xx-xm)/xx;
                Bp = sqrt(Br^2+Bz^2);
                divb=-(xm^2*ym + xx^2*ym - xm^2*yy - xx^2*yy + 3*ym*yy^2 - 3*ym^2*yy + ym^3 - yy^3 - 2*xm*xx*ym + 2*xm*xx*yy)/(xx^4*((xm - xx)^2/xx^2 + (ym - yy)^2/xx^2 + 1)^(3/2));
            case {3,4,5,6,7,9}
                Bt = 1;
                Br = 0;
                Bz = 1;
                Bp = sqrt(Br^2+Bz^2);
                divb = 0;
            case {8,10}
                Br =  (yy-yc);
                Bz = -(xx-xc);
                Bt = 1;
                Bp = sqrt(Br^2+Bz^2);
                divb = 0;
            case 100
                x0 = xc;
                y0 = yc;
                Bt = 1;
                Br = 2*pi/theta*(yy-y0);
                Bz = 2*pi/theta*(-xx+x0);
                Bp = sqrt(Br.^2+Bz.^2);
                divb = 0;
            case 1000
                Bt = 1;
                Br =  (yy-ym)/xx;
                Bz = -(xx-xm)/xx;
                divb=-(xm^2*ym + xx^2*ym - xm^2*yy - xx^2*yy + 3*ym*yy^2 - ...
                    3*ym^2*yy + ym^3 - yy^3 - 2*xm*xx*ym + 2*xm*xx*yy)/(xx^...
                    4*((xm - xx)^2/xx^2 + (ym - yy)^2/xx^2 + 1)^(3/2));
                Bp = sqrt(Br^2+Bz^2);
                if ~axisym
                    error('Only axisym')
                end                
            case 1001
                Bt = 1;
                if axisym
                    Br = 1/xx;
                    Bz = 0;
                    Bp = sqrt(Br^2+Bz^2);
                    divb = 1/(xx^4*(1/xx^2 + 1)^(3/2));
                else
                    error('Only axisym')
                end
            case 1002
                Bt = 1;
                Br = 0;
                Bz = 0;
                Bp = sqrt(Br^2+Bz^2);
                divb = 0;
                
            case 1003
                Bt = 0;
                if axisym
                    Br = 1/xx;
                    Bz = 0;
                    Bp = sqrt(Br^2+Bz^2);
                    divb = 1/(xx^4*(1/xx^2)^(3/2));
                else
                    error('Only axisym')
                end
            case 2000
                Bt = 1;
                Br =  (yy-ym)/xx;
                Bz = -(xx-xm)/xx;
                Bp = sqrt(Br^2+Bz^2);
                divb=-(xm^2*ym + xx^2*ym - xm^2*yy - xx^2*yy + 3*ym*yy^2 - 3*ym^2*yy + ym^3 - yy^3 - 2*xm*xx*ym + 2*xm*xx*yy)/(xx^4*((xm - xx)^2/xx^2 + (ym - yy)^2/xx^2 + 1)^(3/2));
                

        end
        B = sqrt(Bp^2+Bt^2);
        br(i,j) = Br/B;
        bz(i,j) = Bz/B;
        bt(i,j) = Bt/B;
        db(i,j) = divb;
    end
end
b(:,1) = br(:);
b(:,2) = bz(:);
b(:,3) = bt(:);
db = db(:);




