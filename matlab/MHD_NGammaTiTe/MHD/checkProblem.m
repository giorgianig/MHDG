% check if the problem is temporal
if strcmp(problem,'NS_time') || strcmp(problem,'NS_pseudotime')
    % do nothing
elseif strcmp(problem,'NS_steady') || strcmp(problem,'S_steady')
    % set 1 time step and eliminate mass matrix
    nStep = 1;
    M = zeros(size(M));
end

% check if the problem is implicit (NS) or explicit (Stokes)
if strcmp(problem,'NS_time') || strcmp(problem,'NS_steady') || strcmp(problem,'NS_pseudotime')
    % do nothing
elseif strcmp(problem,'S_time') || strcmp(problem,'S_steady')
    % set 1 time step and eliminate mass matrix
    maxIter = 1;
end
