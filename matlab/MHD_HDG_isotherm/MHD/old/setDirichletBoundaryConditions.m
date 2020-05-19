function u = setDirichletBoundaryConditions(X)

neq = 3;
x = X(:,1);
y = X(:,2);
u = zeros(numel(x),neq);

% airfoil or Von Karman Street
% inc = 2;
% ind_outer = x.^2 + y.^2 >2;
% u(ind_outer,1) = 1 * cos(inc*pi/180) ;
% u(ind_outer,2) = 1 * sin(inc*pi/180);

% cavity flow
% OLD
% u(y>0.9,1) = (y(y>0.9)-0.9 )/0.1; 
% NEW
% u(all([x<=0.1,abs(y-1)<1e-6],2),1) = 10*x(all([x<=0.1,abs(y-1)<1e-6],2));
% u(all([x>0.1,x<0.9,abs(y-1)<1e-6],2),1) = 1;
% u(all([x>=0.9,abs(y-1)<1e-6],2),1) = 10-10*x(all([x>=0.9,abs(y-1)<1e-6],2));

% cavity flow
% aux = 2/3;
% u(y>aux,1) = (y(y>aux)-aux )/(1-aux);

% cavity flow
% u(all([x<=0.1,abs(y-1)<1e-6],2),1) = 10*x(all([x<=0.1,abs(y-1)<1e-6],2));
% u(all([x>0.1,x<0.9,abs(y-1)<1e-6],2),1) = 1;
% u(all([x>=0.9,abs(y-1)<1e-6],2),1) = 10-10*x(all([x>=0.9,abs(y-1)<1e-6],2));

%constant orizontal velocity
u(:,1) = 1;
u(:,2) = 0.1;
% u(:,1) = 1;

% turbolence breakdown
% u(:,1) = sin(x*pi).*cos(y*pi);
% u(:,2) = -cos(x*pi).*sin(y*pi);

%% Navier Stokes analytical solutions

% Burger sheet
% x = x-0.5;
% y = y-0.5;
% u(:,1) = 4*x;
% u(:,2) = -4*y;

% Wang 1990a
% x = x - 0.5;
% lambda = 1;
% a = 1;
% b = 1;
% u(:,1) = 2*a*y-b*lambda*exp(-lambda*y).*cos(lambda*x);
% u(:,2) = b*exp(-lambda*y).*sin(lambda*x)*lambda;

% Oblique impingement of two jets
% a = 1;
% b = 1;
% u(:,1) = 2*a*y+b*x;
% u(:,2) = -y*b;

% Tsien 1943a
% x = x+0.5;
% y = y+0.5;
% a = 0.1;
% b = 0.1;
% c = 0.1;
% u(:,1) = a+2*b*y+c./(x.*(1+y.^2./x.^2));
% u(:,2) = c*y./(x.^2.*(1+y.^2./x.^2));

% Tsien 1943b
% x = x-0.5;
% y = y+0.5;
% a = 0.1;
% b = 0.1;
% c = 0.1;
% u(:,1) = a+2*b*y+2*c*y./(x.^2+y.^2);
% u(:,2) = -2*c*x./(x.^2+y.^2);


