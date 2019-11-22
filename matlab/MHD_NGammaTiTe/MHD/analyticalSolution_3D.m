function [u,ux,uy,ut] = analyticalSolution_3D(X,t)

global testcase Mesh eps axisym


N2d = size(X,1);
N1d = size(t,1);

% allocation
u = zeros(N2d*N1d,4);
ux = u;
uy = u;
ut = u;

x = X(:,1);
y = X(:,2);

xmax = Mesh.maxx;
xmin = Mesh.minx;
ymax = Mesh.maxy;
ymin = Mesh.miny;
xm = 0.5*(xmax+xmin);
ym = 0.5*(ymax+ymin);

a = testcase.wavex;
a = 2*pi*a;
b = testcase.wavey;
b = 2*pi*b;
for j=1:N1d
    for  i=1:N2d
        xx = x(i);
        yy = y(i);
%         tt = t(j);        
%         xr = (x-xmin)/(xmax-xmin);
%         yr = (y-ymin)/(ymax-ymin);
%         xx = xr(i);
%         yy = yr(i);
        tt = t(j);
        r = sqrt((xx-xm)^2+(yy-ym)^2);
        ind = (j-1)*N2d+i;
        switch testcase.n
            
            case 1
                
                u(ind,1) = 2+sin(a*xx).*sin(b*yy);
                u(ind,2) = cos(a*xx).*cos(b*yy);
                u(ind,3) = 20+cos(a*xx).*sin(b*yy);
                u(ind,4) = 10-sin(a*xx).*cos(b*yy);
                ux(ind,1) = a*cos(a*xx).*sin(b*yy);
                ux(ind,2) = -a*sin(a*xx).*cos(b*yy);
                ux(ind,3) = -a*sin(a*xx).*sin(b*yy);
                ux(ind,4) = -a*cos(a*xx).*cos(b*yy);
                uy(ind,1) =  b*sin(a*xx).*cos(b*yy);
                uy(ind,2) = -b*cos(a*xx).*sin(b*yy);
                uy(ind,3) =  b*cos(a*xx).*cos(b*yy);
                uy(ind,4) =  b*sin(a*xx).*sin(b*yy);
            case 2
                u(ind,1) = 2+sin(a*xx).*sin(b*yy);
                u(ind,2) = cos(a*xx).*cos(b*yy);
                u(ind,3) = 20+cos(a*xx).*sin(b*yy);
                u(ind,4) = 10-sin(a*xx).*cos(b*yy);
                ux(ind,1) = a*cos(a*xx).*sin(b*yy);
                ux(ind,2) = -a*sin(a*xx).*cos(b*yy);
                ux(ind,3) = -a*sin(a*xx).*sin(b*yy);
                ux(ind,4) = -a*cos(a*xx).*cos(b*yy);
                uy(ind,1) =  b*sin(a*xx).*cos(b*yy);
                uy(ind,2) = -b*cos(a*xx).*sin(b*yy);
                uy(ind,3) =  b*cos(a*xx).*cos(b*yy);
                uy(ind,4) =  b*sin(a*xx).*sin(b*yy);
            case 60
                u(ind,1) = 1;
                u(ind,2) = 0;
                u(ind,3) = 18;
                u(ind,4) = 18;
                
        end
    end
end

% convert to conservative variables
u2x = ( ux(:,1).*u(:,2)+u(:,1).*ux(:,2) );
u2y = ( uy(:,1).*u(:,2)+u(:,1).*uy(:,2) );
u2t = ( ut(:,1).*u(:,2)+u(:,1).*ut(:,2) );
u3x = ( ux(:,1).*u(:,3)+u(:,1).*ux(:,3) );
u3y = ( uy(:,1).*u(:,3)+u(:,1).*uy(:,3) );
u3t = ( ut(:,1).*u(:,3)+u(:,1).*ut(:,3) );
u4x = ( ux(:,1).*u(:,4)+u(:,1).*ux(:,4) );
u4y = ( uy(:,1).*u(:,4)+u(:,1).*uy(:,4) );
u4t = ( ut(:,1).*u(:,4)+u(:,1).*ut(:,4) );
ux(:,2) = u2x; 
uy(:,2) = u2y;
ut(:,2) = u2t;
ux(:,3) = u3x; 
uy(:,3) = u3y;
ut(:,3) = u3t;
ux(:,4) = u4x; 
uy(:,4) = u4y;
ut(:,4) = u4t;
u(:,2) = u(:,1).*u(:,2);
u(:,3) = u(:,1).*u(:,3);
u(:,4) = u(:,1).*u(:,4);

