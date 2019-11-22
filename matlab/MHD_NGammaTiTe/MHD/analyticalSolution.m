function [u,ux,uy] = analyticalSolution(X)

global testcase Mesh eps axisym neq


% u(1) = density;
% u(2) = parallel velocity;

% allocation
u = zeros(size(X,1),neq);
ux = u;
uy = u;

x = X(:,1);
y = X(:,2);

xmax = Mesh.maxx;
xmin = Mesh.minx;
ymax = Mesh.maxy;
ymin = Mesh.miny;


xmax = 1;
xmin = 0;
ymax = 1;
ymin = 0;
switch testcase.n
    
    case 1   
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;    
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
        u(:,3) = 20+cos(a*xr).*sin(b*yr);
        u(:,4) = 10-sin(a*xr).*cos(b*yr);
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        ux(:,3) = -a/(xmax-xmin)*sin(a*xr).*sin(b*yr);
        ux(:,4) = -a/(xmax-xmin)*cos(a*xr).*cos(b*yr);
        uy(:,1) =  b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);
        uy(:,3) =  b/(ymax-ymin)*cos(a*xr).*cos(b*yr);
        uy(:,4) =  b/(ymax-ymin)*sin(a*xr).*sin(b*yr);
        
      case 2 % 
       %% bx = 1/30; by = 1/30; n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
        if axisym,error('This is NOT an axisymmetric case!'),end
        u(:,1) = 1+x.^2+y.^2;
        u(:,2) = x.^2+y.^2; 
        u(:,3) = 10+x.^2+y.^2;
        u(:,4) = 10+x.^2-y.^2;
        ux(:,1) = 2*x;
        ux(:,2) = 2*x;
        ux(:,3) = 2*x;
        ux(:,4) = 2*x;
        uy(:,1) = 2*y;
        uy(:,2) = 2*y;
        uy(:,3) = 2*y;
        uy(:,4) = -2*y;

        
      case 3 % 
       %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n = 1+x^2+y^2, u = x^2+y^2; Ei = 10+x^2+y^2; Ee = 10+x^2-y^2;
        if axisym,error('This is NOT an axisymmetric case!'),end
        u(:,1) = 1+x.^2+y.^2;
        u(:,2) = x.^2+y.^2; 
        u(:,3) = 10+x.^2+y.^2;
        u(:,4) = 10+x.^2-y.^2;
        ux(:,1) = 2*x;
        ux(:,2) = 2*x;
        ux(:,3) = 2*x;
        ux(:,4) = 2*x;
        uy(:,1) = 2*y;
        uy(:,2) = 2*y;
        uy(:,3) = 2*y;
        uy(:,4) = -2*y;
        

      case 4 % 
       %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
        u(:,3) = 20+cos(a*xr).*sin(b*yr);
        u(:,4) = 20+sin(a*xr).*cos(b*yr);
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        ux(:,3) = -a/(xmax-xmin)*sin(a*xr).*sin(b*yr);
        ux(:,4) =  a/(xmax-xmin)*cos(a*xr).*cos(b*yr);
        uy(:,1) =  b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);
        uy(:,3) =  b/(ymax-ymin)*cos(a*xr).*cos(b*yr);
        uy(:,4) = -b/(ymax-ymin)*sin(a*xr).*sin(b*yr);
  
        
      case 5 % test for Bohm bc
       %% bx = 1/30 ; by = 0; n=1.1-sin(xr);u=sin(xr); Ei = 2-sin(xr); Ee = 1.4-sin(xr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        a = pi/2;
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        u(:,1) = 1.1-sin(a*xr);
        u(:,2) = sin(a*xr);
        u(:,3) = 2-sin(a*xr);
        u(:,4) = 1.4-sin(a*xr);
        ux(:,1) = -a/(xmax-xmin)*cos(a*xr);
        ux(:,2) =  a/(xmax-xmin)*cos(a*xr);
        ux(:,3) = -a/(xmax-xmin)*cos(a*xr);
        ux(:,4) = -a/(xmax-xmin)*cos(a*xr);
        uy(:,1) = 0;
        uy(:,2) = 0;
        uy(:,3) = 0;
        uy(:,4) = 0;

        
      case 6 % same as case 5 but to be used with Dirichlet BC
       %% bx = 1/30 ; by = 0; n=1.1-sin(xr);u=sin(xr); Ei = 2-sin(xr); Ee = 1.4-sin(xr);
        if axisym,error('This is NOT an axisymmetric case!'),end
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        u(:,1) = 1.1-sin(xr);
        u(:,2) = sin(xr);
        u(:,3) = 2-sin(xr);
        u(:,4) = 1.4-sin(xr);
        ux(:,1) = -1/(xmax-xmin)*cos(xr);
        ux(:,2) =  1/(xmax-xmin)*cos(xr);
        ux(:,3) = -1/(xmax-xmin)*cos(xr);
        ux(:,4) = -1/(xmax-xmin)*cos(xr);
        uy(:,1) = 0;
        uy(:,2) = 0;
        uy(:,3) = 0;
        uy(:,4) = 0;

    case 25 % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 20+sin(a*xr)*cos(b*yr);
        if ~axisym,error('This is an axisymmetric case!'),end
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
        u(:,3) = 20+cos(a*xr).*sin(b*yr);
        u(:,4) = 20+sin(a*xr).*cos(b*yr);
        ux(:,1) = a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        ux(:,3) = -a/(xmax-xmin)*sin(a*xr).*sin(b*yr);
        ux(:,4) =  a/(xmax-xmin)*cos(a*xr).*cos(b*yr);
        uy(:,1) =  b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);
        uy(:,3) =  b/(ymax-ymin)*cos(a*xr).*cos(b*yr);
        uy(:,4) = -b/(ymax-ymin)*sin(a*xr).*sin(b*yr);

        
    case 26 % Axisimmetric case!
        %% bx = 1/30(x-y.^2+2); by = 1/30(x.*y+y); n=2+sin(a*xr)*sin(b*yr);u=cos(a*xr)*cos(b*yr); 
        %% Ei = 20+cos(a*xr)*sin(b*yr); Ee = 10-sin(a*xr)*cos(b*yr);
        if ~axisym,error('This is an axisymmetric case!'),end
        xr = (x-xmin)/(xmax-xmin);
        yr = (y-ymin)/(ymax-ymin);
        a = testcase.wavex; 
        b = testcase.wavey;        
        a = a*2*pi;
        b = b*2*pi;        
        u(:,1) = 2+sin(a*xr).*sin(b*yr);
        u(:,2) = cos(a*xr).*cos(b*yr);
        u(:,3) = 20+cos(a*xr).*sin(b*yr);
        u(:,4) = 10-sin(a*xr).*cos(b*yr);
        
        ux(:,1) =  a/(xmax-xmin)*cos(a*xr).*sin(b*yr);
        ux(:,2) = -a/(xmax-xmin)*sin(a*xr).*cos(b*yr);
        ux(:,3) = -a/(xmax-xmin)*sin(a*xr).*sin(b*yr);
        ux(:,4) = -a/(xmax-xmin)*cos(a*xr).*cos(b*yr);
        
        uy(:,1) =  b/(ymax-ymin)*sin(a*xr).*cos(b*yr);
        uy(:,2) = -b/(ymax-ymin)*cos(a*xr).*sin(b*yr);
        uy(:,3) =  b/(ymax-ymin)*cos(a*xr).*cos(b*yr);
        uy(:,4) =  b/(ymax-ymin)*sin(a*xr).*sin(b*yr);
        
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
    case {50,51,52,53,54,55}
        % WEST
        r = sqrt(x.^2+y.^2);
        u(:,1) = 1;
        u(:,2) = 0;
        u(:,3) = 36;
        u(:,4) = 36;

%% Circular case with limiter
    case {60,69}
        u(:,1) = 1;
        u(:,2) = 0;
        u(:,3) = 18;
        u(:,4) = 18;
    case 61
        u(:,1) = 1;
        u(:,2) = 0;
        u(:,3) = 36;
        u(:,4) = 36;        
    case 65
        u(:,1) = 1;
        u(:,2) = 0;   
        r = sqrt( (x*Mesh.lscale-3.4).^2+(y*Mesh.lscale-0.75).^2)/Mesh.lscale;
        r0 = 0.05/Mesh.lscale;
        ind = r<r0;
        u(ind,2) = 1;
        
    case 100
        %% 
        u(:,1) = 2+sin(a*x)*sin(a*y);
        u(:,2) = cos(a*x)*cos(a*y);
        u(:,3) = 20+cos(a*x)*sin(a*y);
        u(:,4) = 10-sin(a*x)*cos(a*y);     
end

% convert to conservative variables
u2x = ( ux(:,1).*u(:,2)+u(:,1).*ux(:,2) );
u2y = ( uy(:,1).*u(:,2)+u(:,1).*uy(:,2) );
u3x = ( ux(:,1).*u(:,3)+u(:,1).*ux(:,3) );
u3y = ( uy(:,1).*u(:,3)+u(:,1).*uy(:,3) );
u4x = ( ux(:,1).*u(:,4)+u(:,1).*ux(:,4) );
u4y = ( uy(:,1).*u(:,4)+u(:,1).*uy(:,4) );
ux(:,2) = u2x; 
uy(:,2) = u2y;
ux(:,3) = u3x; 
uy(:,3) = u3y;
ux(:,4) = u4x; 
uy(:,4) = u4y;
u(:,2) = u(:,1).*u(:,2);
u(:,3) = u(:,1).*u(:,3);
u(:,4) = u(:,1).*u(:,4);