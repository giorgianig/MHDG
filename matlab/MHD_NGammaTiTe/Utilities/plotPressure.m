function plotPressure(X,T,p)

nElems = size(T,1);
Np = size(T,2);

% initialize figure
figure(2); clf;
view(3); grid on;
set(gca, 'FontSize', 12,...
    'XTick', [0:0.25:1], 'YTick', [0:0.25:1]);

% loop in elements
for ielem = 1:nElems

    % index
    ind = (ielem-1)*Np + (1:Np);
    Te = T(ielem,:);

    % plot
    tri = delaunay(X(Te,1),X(Te,2));
    trisurf(tri,X(Te,1),X(Te,2),p(ind));
    hold on
end