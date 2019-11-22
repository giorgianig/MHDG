function Dirichlet_boundaries = setDirichletFaces(boundaryNames)
global  testcase

Dirichlet_boundaries = false(numel(boundaryNames),1);

% loop in boundaries
for iboundary = 1:numel(boundaryNames)
    if testcase.n==69 || testcase.n==26
        Dirichlet_boundaries(iboundary) = true;
    else
        iname = boundaryNames{iboundary}(4:end);
        if testcase.n==5
            if any(strcmp(iname,{'Diriclet','UP','DOWN','LEFT'}))
                Dirichlet_boundaries(iboundary) = true;
            end
        else
            if any(strcmp(iname,{'Diriclet','UP','DOWN','LEFT','RIGHT'}))
                Dirichlet_boundaries(iboundary) = true;
            end
        end
    end
end

