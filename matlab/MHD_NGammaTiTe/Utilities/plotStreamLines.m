function plotStreamLines(X,T,u,np,sx,sy)
% u is the continue solution!

% velocity components
u_x = u(1:2:end-1);
u_y = u(2:2:end);

% plotting grid
x = linspace(min(X(:,1)),max(X(:,1)),np);
y = linspace(min(X(:,2)),max(X(:,2)),np);
[X_plot,Y_plot] = meshgrid(x,y);

% Velocity in the plotting mesh
Ux_plot = griddata(X(:,1),X(:,2),u_x,X_plot,Y_plot,'cubic');
Uy_plot = griddata(X(:,1),X(:,2),u_y,X_plot,Y_plot,'cubic');

% plot streamlines
figure
plotMesh(X,T)
streamline(X_plot,Y_plot,Ux_plot,Uy_plot,sx,sy);
axis equal
