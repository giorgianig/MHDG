function res = filterSolutionElement(u,refEl)


coordRef = refEl.NodesCoord;
np = size(coordRef,1);
nDeg = refEl.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = V\eye(np);

% Convert solution into modal expansion
um = invV*u;

npm1 = nDeg*(nDeg+1)/2;
% umho = zeros(size(um));
% umho(npm1+1:np,:) = um(npm1+1:np,:);
ulo = zeros(size(um));
ulo(1) = um(1);
res = V*ulo;
