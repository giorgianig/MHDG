function eps = findCoeffShockCaptur(u,X,T,refEl,k,neq,thresh)

ne = size(T,1); % number of elements
np = size(T,2); % number of nodes x element
u = transpose(reshape(u,2,numel(u)/neq));

% use rho for shock detection
% ud = u(:,1);
% use gamma for shock detection
% ud = u(:,2);
% use u_par for shock detection
ud = u(:,2)./u(:,1);
ud = reshape(ud,np,ne);

% Vandermonde matrix
coordRef = refEl.NodesCoord;
nDeg = refEl.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = V\eye(np);

% Convert solution into modal expansion
um = invV*ud;

% Solution with only the ho mode
npm1 = nDeg*(nDeg+1)/2;
umho = zeros(size(um));
umho(npm1+1:np,:) = um(npm1+1:np,:);

% Shock detector
se = log10( sum(umho.^2,1)./sum(um.^2,1))';
indth = sum(um.^2,1)<thresh;
se(indth)= -inf;

% coefficients
s0 = log10(1/nDeg^4);
eps0 = 1e-2*findElementSize(X,T)/nDeg;

% set epsilon
ind1 = all([se>=s0-k,se<=s0+k],2);
ind2 = se>s0+k;
eps = zeros(ne,1);
eps(ind1) = eps0(ind1)/2.*(1+sin(pi/(2*k)*(se(ind1)-s0)));
eps(ind2) = eps0(ind2);


