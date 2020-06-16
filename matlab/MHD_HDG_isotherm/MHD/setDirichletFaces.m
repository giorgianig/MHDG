function Dirichlet_boundaries = setDirichletFaces(boundaryNames)

Dirichlet_boundaries = false(numel(boundaryNames),1);

% loop in boundaries
for iboundary = 1:numel(boundaryNames)
    iname = boundaryNames{iboundary}(4:end);
        if any(strcmp(iname,{'Diriclet'}))
            Dirichlet_boundaries(iboundary) = true;
        end
end

