function [b,db,dvx,dvy] = defineMagneticField(X)

global testcase Magnetic Mesh axisym driftvel
% allocation
b = zeros(size(X));
db = b; % dbx_dx, dby_dy

x = X(:,1);
y = X(:,2);

dvx = b;
dvy = b;

xc = testcase.xc;
yc = testcase.yc;
xmax = max(Mesh.X(:,1));
xmin = min(Mesh.X(:,1));
ymax = max(Mesh.X(:,2));
ymin = min(Mesh.X(:,2));
xm = 0.5*(xmax+xmin);
ym = 0.5*(ymax+ymin);
xm = xc;
ym = yc;
%% bx = 0.1*Bx/|B| -- by = 0.1*By/|B|
switch testcase.n
    case 1
        for i=1:size(x)
            xx = x(i);
            yy = y(i);
            Br =  (yy-ym);
            Bz = -(xx-xm);
            Bt = 1;
            B = sqrt(Br^2+Bz^2+Bt^2);
            b(i,1) = Br/B;
            b(i,2) = Bz/B;
        end
        % db = 0;  
    case 2
       %% bx = 1/30; by = 1/30; n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
        b(:,1) =1/30;
        b(:,2) = 1/30;
        % db = 0;  
    case 3
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
        if axisym,error('This is NOT an axisymmetric case!'),end
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30;
        db(:,2) = 1/30*(x+1);  
    case 4
       %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30;
        db(:,2) = 1/30*(x+1);          
      case 5 % test for Bohm bc
       %% bx = 1/30 ; by = 0; n=1.1-sin(xr);u=sin(xr); Ei = 2-sin(xr); Ee = 1.4-sin(xr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        b(:,1) = 1/30;
        b(:,2) = 0;
      case 6 % test for Bohm bc
       %% bx = 1/30 ; by = 0; n=1.1-sin(xr);u=sin(xr); Ei = 2-sin(xr); Ee = 1.4-sin(xr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        b(:,1) = 1/30;
        b(:,2) = 0;        
%     case 24
%        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
%         if ~axisym,error('This is an axisymmetric case!'),end
%         b(:,1) = 1/30*(x-y.^2+2);
%         b(:,2) = 1/30*(x.*y+y);
%         db(:,1) = 1/30+((1/30)*x-(1/30)*y.^2+1/15)./x;
%         db(:,2) = 1/30*(x+1);          
    case 25 % Axisimmetric case!
      %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        if ~axisym,error('This is an axisymmetric case!'),end
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30+((1/30)*x-(1/30)*y.^2+1/15)./x;
        db(:,2) = 1/30*(x+1);          
    case 26 % Axisimmetric case!
      %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20-sin(a*xr)*cos(b*yr);
        if ~axisym,error('This is an axisymmetric case!'),end
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30+((1/30)*x-(1/30)*y.^2+1/15)./x;
        db(:,2) = 1/30*(x+1);             
    case 30
        %% Circular field centered in [0, 0], n: density spot, u=0
        b(:,1) = 0.1*y./sqrt(y.^2+x.^2);
        b(:,2) = 0.1*-x./sqrt(y.^2+x.^2);
        db(:,1) = 0.1* ( -y.*x./(y.^2+x.^2).^(3/2) );
        db(:,2) = 0.1* (  y.*x./(y.^2+x.^2).^(3/2) ) ;
        indzero = (x.^2+y.^2)<1e-6;
        b(indzero,:) = 0;
        db(indzero,:) = 0;
    case 31
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
        b(:,1) = 0.1*y./sqrt(y.^2+x.^2);
        b(:,2) = 0.1*-x./sqrt(y.^2+x.^2);
        db(:,1) = 0.1* ( -y.*x./(y.^2+x.^2).^(3/2) );
        db(:,2) = 0.1* (  y.*x./(y.^2+x.^2).^(3/2) ) ;
        indzero = (x.^2+y.^2)<1e-6;
        b(indzero,:) = 0;
        db(indzero,:) = 0;        
    case 32
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
        b(:,1) = 0.1*y./sqrt(y.^2+x.^2);
        b(:,2) = 0.1*-x./sqrt(y.^2+x.^2);
        db(:,1) = 0.1* ( -y.*x./(y.^2+x.^2).^(3/2) );
        db(:,2) = 0.1* (  y.*x./(y.^2+x.^2).^(3/2) ) ;
        indzero = (x.^2+y.^2)<1e-6;
        b(indzero,:) = 0;
        db(indzero,:) = 0;        
        
%% West        
    case 50
        % WEST
        error('You shouldnt be here!')
    case 51
        % WEST
        error('You shouldnt be here!')
    case 52
        % WEST
        error('You shouldnt be here!')        
    case 53
        % WEST
        error('You shouldnt be here!')   
    case 54
        % WEST
        error('You shouldnt be here!')   
    case 55
        % WEST
        error('You shouldnt be here!')           
 %% Circular case with limiter
    case {60,69} 
        
        x = x*Mesh.lscale;
        y = y*Mesh.lscale;        
        R0 = 0.5*(Mesh.maxx+Mesh.minx)*Mesh.lscale;        
        q = 3.5;
        r = sqrt((x-R0).^2+y.^2);
        b(:,1) = -y./sqrt(R0^2*q^2+(1-q^2)*r.^2);
        b(:,2) = (x-R0)./sqrt(R0^2*q^2+(1-q^2)*r.^2);
        if axisym
            db(:,1) = -y./x./sqrt(R0^2*q^2+(1-q^2)*r.^2)*Mesh.lscale;
        end
        
        if (nargout>2 && driftvel)
            
%             Bm = sqrt((y.^2+x.^2-2*x*R0+R0^2-q^2*x.^2+2*q^2*x*R0-q^2*y.^2)./(-x.^2+2*x*R0-y.^2))./x;
            
%             relx = (1/2)*((2*y-2*q^2*y)./(-x.^2+2*x*R0-y.^2)+(2*(y.^2+x.^2-2*x*R0+R0^2-q^2*x.^2+...
%                 2*q^2*x*R0-q^2*y.^2)).*y./(-x.^2+2*x*R0-y.^2).^2)./(x.*sqrt((y.^2+x.^2-2*x*R0+R0^2-...
%                 q^2*x.^2+2*q^2*x*R0-q^2*y.^2)./(-x.^2+2*x*R0-y.^2)));
%             
%             rely = -sqrt((y.^2+x.^2-2*x*R0+R0^2-q^2*x.^2+2*q^2*x*R0-q^2*y.^2)./(-x.^2+2*x*R0-y.^2))./x.^2+...
%                 (1/2)*((2*x-2*R0-2*q^2*x+2*q^2*R0)./(-x.^2+2*x*R0-y.^2)-(y.^2+x.^2-2*x*R0+R0^2-...
%                 q^2*x.^2+2*q^2*x*R0-q^2*y.^2).*(-2*x+2*R0)./(-x.^2+2*x*R0-y.^2).^2)./(x.*sqrt((y.^2+x.^2-...
%                 2*x*R0+R0^2-q^2*x.^2+2*q^2*x*R0-q^2*y.^2)./(-x.^2+2*x*R0-y.^2)));
%             
%             dvx = q./x.*relx./Bm.^3;
%             dvy =-q./x.*rely./Bm.^3;
            
            
            
            Bm = R0./x;
            rely = -R0./x.^2;
            relx = 0;
            dvx = -Bm.*relx./Bm.^3*Mesh.lscale;
            dvy = Bm.*rely./Bm.^3*Mesh.lscale;            
            
        end
        
    case 100
        %% Bx = 1, By = , n = 1, u = 0, E = 1
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30;
        db(:,2) = 1/30+(1/30)*x;
              
        
        
end






