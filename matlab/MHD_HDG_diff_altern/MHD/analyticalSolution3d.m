function [u,ux,uy,ut] = analyticalSolution3d(X,t)

global testcase Mesh eps axisym


N2d = size(X,1);
N1d = size(t,1);

% allocation
u = zeros(N2d*N1d,2);
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
        tt = t(j);
        r = sqrt((xx-xm)^2+(yy-ym)^2);
        ind = (j-1)*N2d+i;
        switch testcase.n
            
            case 1
                %% Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 2+sin(a* xx)*sin(a* yy);
                u(ind,2) = cos(a* xx)*cos(a* yy);
                ux(ind,1) = a*cos(a* xx)*sin(a* yy);
                ux(ind,2) = -a*sin(a* xx)*cos(a* yy);
                uy(ind,1) = a*sin(a* xx)*cos(a* yy);
                uy(ind,2) = -a*cos(a* xx)*sin(a* yy);
            case 2
                if ~axisym 
                    error("This is an axisymmatric case!")
                end
                %% ! Axisimmetric case with div(b)~=0,circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 2+sin(a* xx)*sin(a* yy);
                u(ind,2) = cos(a* xx)*cos(a* yy);
                ux(ind,1) = a*cos(a* xx)*sin(a* yy);
                ux(ind,2) = -a*sin(a* xx)*cos(a* yy);
                uy(ind,1) = a*sin(a* xx)*cos(a* yy);
                uy(ind,2) = -a*cos(a* xx)*sin(a* yy);                
            case 3
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1+xx;
                u(ind,2) = yy;
                ux(ind,1) = 1;
                ux(ind,2) = 0;
                uy(ind,1) = 0;
                uy(ind,2) = 1;   
            case 4
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1;
                u(ind,2) = xx;
                ux(ind,1) = 0;
                ux(ind,2) = 1;
                uy(ind,1) = 0;
                uy(ind,2) = 0;            
            case 5
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1;
                u(ind,2) = 1;
                ux(ind,1) = 0;
                ux(ind,2) = 0;
                uy(ind,1) = 0;
                uy(ind,2) = 0;      
            case 6
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1+xx;
                u(ind,2) = xx^3;
                ux(ind,1) = 1;
                ux(ind,2) = 3*xx^2;
                uy(ind,1) = 0;
                uy(ind,2) = 0;         
            case 7
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1+xx+yy;
                u(ind,2) = xx^2+yy^2;
                ux(ind,1) = 1;
                ux(ind,2) = 2*xx;
                uy(ind,1) = 1;
                uy(ind,2) = 2*yy;         
            case 8
                %% ! Cartesian case, , circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 1+xx;
                u(ind,2) = 1+yy;
                ux(ind,1) = 1;
                ux(ind,2) = 0;
                uy(ind,1) = 0;
                uy(ind,2) = 1;   
            case 9
                %% ! Cartesian case, Bx = 0, Bz = 1, Bt = 1
                u(ind,1) = 1+xx;
                u(ind,2) = 1+yy;
                ux(ind,1) = 1;
                ux(ind,2) = 0;
                uy(ind,1) = 0;
                uy(ind,2) = 1;     
            case 10
                %% ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 1;
                u(ind,2) = 1;
                ux(ind,1) = 0;
                ux(ind,2) = 0;
                uy(ind,1) = 0;
                uy(ind,2) = 0;     
            case 100
                %% Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 2+sin(a* xx)*sin(a* yy)*sin(b*tt);
                u(ind,2) = cos(a* xx)*cos(a* yy)*cos(b*tt);                
                ux(ind,1) = a*cos(a* xx)*sin(a* yy)*sin(b*tt);
                ux(ind,2) = -a*sin(a* xx)*cos(a* yy)*cos(b*tt);
                uy(ind,1) = a*sin(a* xx)*cos(a* yy)*sin(b*tt);
                uy(ind,2) = -a*cos(a* xx)*sin(a* yy)*cos(b*tt);
                ut(ind,1) = b*sin(a* xx)*sin(a* yy)*cos(b*tt);
                ut(ind,2) = -b*cos(a* xx)*cos(a* yy)*sin(b*tt);  
            case 1000
                u(ind,1) = 1;
                u(ind,2) = 1;
            case 1001
                u(ind,1) = 1;
                u(ind,2) = 1;
            case 1002
                u(ind,1) = 1;
                u(ind,2) = 1;                
            case 1003
                u(ind,1) = 1;
                u(ind,2) = 1;                  
            case 2000
                if ~axisym 
                    error("This is an axisymmatric case!")
                end
                %% ! Axisimmetric case with div(b)~=0,circular field centered in [xm, ym] in the poloidal plane, Bt = 1
                u(ind,1) = 1;
                u(ind,2) = cos(a* xx)*cos(a* yy);
                ux(ind,1) = 0;
                ux(ind,2) = -a*sin(a* xx)*cos(a* yy);
                uy(ind,1) = 0;
                uy(ind,2) = -a*cos(a* xx)*sin(a* yy);                   
                
        end
    end
end

% convert to conservative variables
u2x = ( ux(:,1).*u(:,2)+u(:,1).*ux(:,2) );
u2y = ( uy(:,1).*u(:,2)+u(:,1).*uy(:,2) );
u2t = ( ut(:,1).*u(:,2)+u(:,1).*ut(:,2) );
ux(:,2) = u2x;
uy(:,2) = u2y;
ut(:,2) = u2t;
u(:,2) = u(:,1).*u(:,2);


