function [u,ux,uy] = analyticalSolution(X)

global testcase Mesh eps axisym


% u(1) = density;
% u(2) = parallel velocity;

% allocation
u = zeros(size(X));
ux = u;
uy = u;

x = X(:,1);
y = X(:,2);

xmax = Mesh.maxx;
xmin = Mesh.minx;
ymax = Mesh.maxy;
ymin = Mesh.miny;

switch testcase.n
    
    case 1
        %% Bx = 1, By = 1, n = 1, u = 1
        u(:,1) = 1;
        u(:,2) = 1;
    case 2
        %% Bx = 1, By = 1, n = x+1, u = 1 
        u(:,1) = x+1;
        u(:,2) = 1;
        ux(:,1) = 1;
    case 3
        %% Bx = 1, By = 1, n = 1, u = x 
        u(:,1) = 1;
        u(:,2) = x;    
        ux(:,2) = 1;
    case 4        
        %% Bx = 1, By = 1, n = x+1, u = y
        u(:,1) = x+1;
        u(:,2) = y;
        ux(:,1) = 1;
        uy(:,2) = 1;
    case 5
        %% Bx = y, By = 1, n = 1, u = 1
        u(:,1) = 1;
        u(:,2) = 1;
    case 6
        %% Circular field centered in [0, 0], n=1, u=1
        u(:,1) = 1;
        u(:,2) = 1;    
     case 7
        %% Circular field centered in [0, 0], n=x+1, u=1
        u(:,1) = x+1;
        u(:,2) = 1; 
        ux(:,1) = 1;
     case 8
        %% Circular field centered in [0, 0], n=x+1, u=y
        u(:,1) = x+1;
        u(:,2) = y; 
        ux(:,1) = 1;
        uy(:,2) = 1;
    case 9
        %% Circular field centered in [-0.5, -0.5], n = 2+sin(wx*x )*sin(wy*y)
        %% u = cos(wx*x)*cos(wy*y)
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);

        a = testcase.wavex; 
        b = testcase.wavey;
        
        a = a*2*pi;
        b = b*2*pi;
        u(:,1) = 2 + sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        uy(:,1) = b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);
    case 10
        %% Bx = 1, By = 1, n = 2+sin(wx*x )*sin(wy*y)
        %% u = cos(wx*x)*cos(wy*y)
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);

        a = testcase.wavex; 
        b = testcase.wavey;
        
        a = a*2*pi;
        b = b*2*pi;
        u(:,1) = 2 + sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);  
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        uy(:,1) = b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);        
       
    case 11
        %% Bx = 1, By = 1, n = 1
        %% u = cos(wx*x)*cos(wy*y)
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);

        a = testcase.wavex; 
        b = testcase.wavey;
        
        a = a*2*pi;
        b = b*2*pi;
        u(:,1) = 1;
        u(:,2) = cos(a*xr).*cos(b*yr);    
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);        
    case 12
        %% Bx = 1, By = 1, n = 1, u = x^2+y^2 
        u(:,1) = 1;
        u(:,2) = x.^2+y.^2;     
        ux(:,2) = 2*x;
        uy(:,2) = 2*y;
    case 13
        %% Bx = 1, By = 1, n = 1
        %% u = cos(wx*x)
        xr = (x-xmin)/(xmax-xmin);
        a = testcase.wavex;                 
        a = a*2*pi;
        u(:,1) = 1;
        u(:,2) = cos(a*xr);            
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr);
    case 14
        %% Bx = 1, By = 1, n = 1+x+y, u = x^2+y^2        
        u(:,1) = 1+x+y;
        u(:,2) = x.^2+y.^2;       
        ux(:,1) = 1;
        ux(:,2) = 2*x;
        uy(:,1) = 1;
        uy(:,2) = 2*y;
    case 15
        %% Bx = 1, By = 1, n = 2+cos(wx*x), u = 1
        xr = (x-xmin)/(xmax-xmin);
        a = testcase.wavex;                 
        a = a*2*pi;
        u(:,1) = 2+cos(a*xr); 
        u(:,2) = 1;   
        ux(:,1) = -a/(xmax-xmin)*sin(a*xr);
    case 16
        %% Bx = 1, By = 1, n = 2+sin(wx*x)*sin(wy*y)
        %% u = 1
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = 1;
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        uy(:,1) = b/(ymax-ymin)*sin(a*xr).*cos(b*yr);

    case 17
        %% Circular field centered in [xc, yc], n = 2+x^2+y^2
        %% u = x^3+y^3
        u(:,1) = 2+x.^2+y.^2;
        u(:,2) = x.^3+y.^3;     
        ux(:,1) = 2*x;
        uy(:,1) = 2*y;        
        ux(:,2) = 3*x.^2;     
        uy(:,2) = 3*y.^2;    
    case 18
        %%  Circular field centered in [xc, yc], n = 1+x+y, u = x^2+y^2        
        u(:,1) = 1+x+y;
        u(:,2) = x.^2+y.^2;       
        ux(:,1) = 1;
        ux(:,2) = 2*x;
        uy(:,1) = 1;
        uy(:,2) = 2*y; 
    case 23
        % bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=1;u=1;
        u(:,1) = 1;
        u(:,2) = 1;       
        
    case 24 
        % bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
    case 25 % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr);
        if ~axisym,error('This is an axisymmetric case!'),end
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
%     case 20
%         %%  Circular field centered in [xc, yc], n = 1+x+y, u = x+y 
%         u(:,1) = 1+x+y;
%         u(:,2) = x+y;       
%         ux(:,1) = 1;
%         ux(:,2) = 1;
%         uy(:,1) = 1;
%         uy(:,2) = 1;
%     case 21
%         %% Bx = 1, By = 1, n = 1+x^2+y^2
%         %% u = x^3+y^3
%         u(:,1) = 1+x.^2+y.^2;
%         u(:,2) = x.^3+y.^3;     
%         ux(:,1) = 2*x;
%         uy(:,1) = 2*y;        
%         ux(:,2) = 3*x.^2;     
%         uy(:,2) = 3*y.^2;   
%     case 22
%         %% Bx = -1, By = 1, n = 1+x^2+y^2
%         %% u = x^3+y^3
%         u(:,1) = 1+x.^2+y.^2;
%         u(:,2) = x.^3+y.^3;     
%         ux(:,1) = 2*x;
%         uy(:,1) = 2*y;        
%         ux(:,2) = 3*x.^2;     
%         uy(:,2) = 3*y.^2;         
%         
    case 30
        %% Circular field centered in [0, 0], n: density spot, u=0
        R = 0.7;
        delta = 0.1;
        T = 5; % degrees
        T = T/180*pi;
        [theta,r] = cart2pol(x,y);
        ind = all([abs(theta)<T r<R+delta r>R-delta],2);
        u(:,1) = 1;
        u(ind,1) = 1.1 ; 
        u(:,2) = 0;   
    case 31
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
        R = 0.7;
        delta = 0.05;
        T = 5; % degrees
        T = T/180*pi;
        m = 4;
        [theta,r] = cart2pol(x,y);
        u(:,1) = 1+ m*exp(- ((r-R).^2+(theta-T).^2)/delta^2);
        u(:,2) = 0;  
    case 32
        %% Circular field centered in [0, 0], n: density spot Gaussian, u=0
        R = 0.7;
        delta = 0.05;
        T = 5; % degrees
        T = T/180*pi;
        m = 1;
        [theta,r] = cart2pol(x,y);
        u(:,1) = eps+ m*exp(- ((r-R).^2+(theta-T).^2)/delta^2);
        u(:,2) = 0;  
        
%% West
    case 50
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;
    case 51
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;   
    case 52
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;     
    case 53
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;        
    case 54
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;
    case 55
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;            
%% Circular case with limiter
    case 60
        u(:,1) = 1;
        u(:,2) = 0;   
    case 65
        u(:,1) = 1;
        u(:,2) = 0;   
        r = sqrt( (x*Mesh.lscale-3.4).^2+(y*Mesh.lscale-0.75).^2)/Mesh.lscale;
        r0 = 0.05/Mesh.lscale;
        ind = r<r0;
        u(ind,2) = 1;
end

% convert to conservative variables
u2x = ( ux(:,1).*u(:,2)+u(:,1).*ux(:,2) );
u2y = ( uy(:,1).*u(:,2)+u(:,1).*uy(:,2) );
ux(:,2) = u2x; 
uy(:,2) = u2y;
u(:,2) = u(:,1).*u(:,2);