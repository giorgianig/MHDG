function plotBoundary(X,boundary)

% plot only the boundary of a mesh
boundaryNames = fieldnames(boundary);
nOfBoundaries = numel(boundaryNames);

for i=1:nOfBoundaries
    name = boundaryNames{i};
    Tb = boundary.(name);
    
    hold on
    for j=1:size(Tb,1)
        Tf = Tb(j,:);
        Xf = X(Tf,:);
        plot(Xf(:,1),Xf(:,2),'k-','LineWidth',1);
    end
    hold off
end
axis equal