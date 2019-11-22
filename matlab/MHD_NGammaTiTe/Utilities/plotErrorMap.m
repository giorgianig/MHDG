function plotErrorMap(X,T,error)

% linearize T
T = T(:,1:3);

%Plot
patch('Faces',T,'Vertices',X,'FaceVertexCData',error,...
    'FaceColor','flat','EdgeAlpha',0);
axis equal
colorbar('location','East');