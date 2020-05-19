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

%% bx = 0.1*Bx/|B| -- by = 0.1*By/|B|
switch testcase.n
    case 1
        %% Bx = 1, By = 1, n = 1, u = 1
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 2
        %% Bx = 1, By = 1, n = x+1, u = 1
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 3
        %% Bx = 1, By = 1, n = 1, u = x
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 4
        %% Bx = 1, By = 1, n = x+1, u = y
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 5
        %% Bx = y, By = 1, n = 1, u = 1
        b(:,1) = 0.1*y./sqrt(y.^2+1);
        b(:,2) = 0.1*1./sqrt(y.^2+1);
        db(:,1) = 0;
        db(:,2) = -y./(y.^2+1).^(3/2)*0.1;
    case 6
        %% Circular field centered in [xc, yc], n=1, u=1      
        for i = 1:numel(x)
            Br =  (y(i)-yc);
            Bz = -(x(i)-xc);
            Bt = 1;
            Bp = sqrt(Br^2+Bz^2);
            B = sqrt(Bp^2+Bt^2);
            b(i,1) = Br/B;
            b(i,2) = Bz/B;
%             b(i,1) = 0.1*(y(i)-yc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
%             b(i,2) = 0.1*(-x(i)+xc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
%             db(i,1) = -(1/2)*0.1*(y(i)-yc)*(2*x(i)-2*xc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);
%             db(i,2) = -(1/2)*0.1*(-x(i)+xc)*(2*y(i)-2*yc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);            
        end
    case 7
        %% Circular field centered in [xc, yc], n=x+1, u=1
        for i = 1:numel(x)
            b(i,1) = 0.1*(y(i)-yc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
            b(i,2) = 0.1*(-x(i)+xc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y(i)-yc)*(2*x(i)-2*xc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x(i)+xc)*(2*y(i)-2*yc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);            
        end
    case 8
        %% Circular field centered in [xc, yc], n=x+1, u=y
        for i = 1:numel(x)
            b(i,1) = 0.1*(y(i)-yc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
            b(i,2) = 0.1*(-x(i)+xc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y(i)-yc)*(2*x(i)-2*xc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x(i)+xc)*(2*y(i)-2*yc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);            
        end
    case 9
        %% Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y)
        %% u = cos(wx*x)*cos(wy*y)
        for i = 1:numel(x)
            Br =  (y(i)-yc);
            Bz = -(x(i)-xc);
            Bt = 1;
            Bp = sqrt(Br^2+Bz^2);
            B = sqrt(Bp^2+Bt^2);
            b(i,1) = Br/B;
            b(i,2) = Bz/B;
            
%             b(i,1) = 0.1*(y(i)-yc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
%             b(i,2) = 0.1*(-x(i)+xc)/sqrt(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2);
%             db(i,1) = -(1/2)*0.1*(y(i)-yc)*(2*x(i)-2*xc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);
%             db(i,2) = -(1/2)*0.1*(-x(i)+xc)*(2*y(i)-2*yc)/(y(i)^2-2*y(i)*yc+yc^2+x(i)^2-2*x(i)*xc+xc^2)^(3/2);            
        end
    case 10
        %% Bx = 1, By = 1, n = 2+sin(wx*x )*sin(wy*y)
        %% u = cos(wx*x)*cos(wy*y)
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
        
    case 11
        %% Bx = 1, By = 1, n = 1
        %% u = cos(wx*x)*cos(wy*y)
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 12
        %% Bx = 1, By = 1, n = 1, u = x.^2+y.^2
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 13
        %% Bx = 1, By = 1, n = 1, u = cos(wx*x)
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 14
        %% Bx = 1, By = 1, n = 1+x+y, u = x.^2+y.^2
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 15
        %% Bx = 1, By = 1, n = 2+cos(wx*x), u = 1
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 16
        %% Bx = 1, By = 1, n = 2+sin(wx*x)*sin(wy*y)
        %% u = 1
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;
        % db = 0;
    case 17
        %% Circular field centered in [xc, yc], n = 2+x^2*y^2
        %% u = x^3*y^3
        for i = 1:numel(x)
            b(i,1) = 0.1*(y-yc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            b(i,2) = 0.1*(-x+xc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y-yc)*(2*x-2*xc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x+xc)*(2*y-2*yc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);            
        end        
    case 18
        %% Circular field centered in [xc, yc],n = 2+x^2+y^2
        %% u = x^3+y^3
        for i = 1:numel(x)
            b(i,1) = 0.1*(y-yc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            b(i,2) = 0.1*(-x+xc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y-yc)*(2*x-2*xc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x+xc)*(2*y-2*yc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);            
        end        
    case 19
        %% Circular field centered in [xc, yc],n = 1+x+y, u = x.^2+y.^2
        for i = 1:numel(x)
            b(i,1) = 0.1*(y-yc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            b(i,2) = 0.1*(-x+xc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y-yc)*(2*x-2*xc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x+xc)*(2*y-2*yc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);            
        end   
    case 20
        %% Circular field centered in [xc, yc],n = 1+x+y, u = x+y
        for i = 1:numel(x)
            b(i,1) = 0.1*(y-yc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            b(i,2) = 0.1*(-x+xc)/sqrt(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2);
            db(i,1) = -(1/2)*0.1*(y-yc)*(2*x-2*xc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);
            db(i,2) = -(1/2)*0.1*(-x+xc)*(2*y-2*yc)/(y^2-2*y*yc+yc^2+x^2-2*x*xc+xc^2)^(3/2);            
        end    
    case 21
        %% Bx = 1, By = 1, n = 1+x^2+y^2, u = x^3+y^3
        b(:,1) = 0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;        
    case 22
        %% Bx = -1, By = 1, n = 1+x^2+y^2, u = x^3+y^3
        b(:,1) = -0.1*sqrt(2)/2;
        b(:,2) = 0.1*sqrt(2)/2;   
    case 23
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=1;u=1;
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30;
        db(:,2) = 1/30*(x+1);
    case 24
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        b(:,1) = 1/30*(x-y.^2+2);
        b(:,2) = 1/30*(x.*y+y);
        db(:,1) = 1/30;
        db(:,2) = 1/30*(x+1);      
    case 25 % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
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
    case 60 
        
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
        
        
        
end






