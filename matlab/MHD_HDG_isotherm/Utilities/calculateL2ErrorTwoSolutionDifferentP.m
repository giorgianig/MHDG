function [error eg domArea] = calculateL2ErrorTwoSolutionDifferentP(X1,T1,u1,refEl1,X2,T2,u2,refEl2)

% calculate the L2 error element by element for two meshes with different p

if refEl1.degree >= refEl2.degree
    T_P = T1;
    X_P = X1;
    u_P = u1;
    referenceElement_P = refEl1;
    T_p = T2;
    X_p = X2;
    u_p = u2;
    referenceElement_p = refEl2;
else
    T_P = T2;
    X_P = X2;
    u_P = u2;
    referenceElement_P = refEl2;
    T_p = T1;
    X_p = X1;
    u_p = u1;
    referenceElement_p = refEl1;
end
nOfElements = size(T_P,1);

% number of nodes per element
nOfNodes_p = size(referenceElement_p.NodesCoord,1); % N of nodes lower mesh
nOfNodes_P = size(referenceElement_P.NodesCoord,1); % N of nodes higher mesh

% reshape velocity
u_p = transpose(reshape(u_p,2,nOfElements*nOfNodes_p));
u_P = transpose(reshape(u_P,2,nOfElements*nOfNodes_P));

% Compute shape functions at Gauss points of the higher mesh
coordGauss = referenceElement_P.IPcoordinates; % Gauss points higher mesh
coordRef = referenceElement_p.NodesCoord; % Nodes on the lower mesh
nOfGauss = size(coordGauss,1);
nDeg = referenceElement_p.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');
shapeFunctions = zeros(nOfGauss,nOfNodes_p);
for ipoint = 1:nOfGauss
    [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(coordGauss(ipoint,:),nDeg);
    N = (invV*[p,dp_dxi,dp_deta])';
    shapeFunctions(ipoint,:,1) = N(1,:);
    shapeFunctions(ipoint,:,2) = N(2,:);
    shapeFunctions(ipoint,:,3) = N(3,:);
end

% allocate vectors
error = zeros(nOfElements,1);
solNorm = zeros(nOfElements,1);
domArea = ones(nOfElements,1);
% Loop in elements
for ielem = elements2compute
    Te_P = T_P(ielem,:);
    Xe_P = X_P(Te_P,:);
    Te_p = T_p(ielem,:);
    Xe_p = X_p(Te_p,:);
    ind_p =  (ielem-1)*nOfNodes_p+1:ielem*nOfNodes_p;
    ue_p = u_p(ind_p,:);
    ind_P =  (ielem-1)*nOfNodes_P+1:ielem*nOfNodes_P;
    ue_P = u_P(ind_P,:);
    [elemError elemSolNorm elemArea] = calculateElementalError...
        (ue_P,ue_p,referenceElement_P,Xe_P,Xe_p,shapeFunctions,...
        coordGauss,coordRef,nDeg);
    domArea(ielem) = elemArea;
    solNorm(ielem) = elemSolNorm;
    error(ielem) = elemError;
end

error = sqrt(error2/sum(solNorm));
eg = sqrt(sum(error2)/sum(solNorm));

function [elemError solNorm elemArea] = calculateElementalError...
    (Ue,ue,referenceElement,Xe_P,Xe_p,n,coordGauss,coordRef,nDeg)

%Information of the reference element
IPw = referenceElement.IPweights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;

xyg_xi = Nxi*Xe_P;
xyg_eta = Neta*Xe_P;
detJ = xyg_xi(:,1).*xyg_eta(:,2) - xyg_xi(:,2).*xyg_eta(:,1);
dvolu = detJ.*IPw;

% check gauss point position
xyg_P = N*Xe_P;
xyg_p = n(:,:,1)*Xe_p;
xyg_err = (xyg_P - xyg_p);
if max(abs(xyg_err(:)))>1e-10
    gxy0 = transpose(coordGauss);
    [gxy,n,convergence] = recomputeGaussPoint(Xe_p,xyg_P,nDeg,coordRef,gxy0,n);
    if ~convergence
        disp('not converging')
    end
end
Ueg = N*Ue;
ueg = n(:,:,1)*ue;
err_gauss = Ueg-ueg;
elemError = err_gauss(:,1)'*(dvolu.*err_gauss(:,1)) + err_gauss(:,2)'*(dvolu.*err_gauss(:,2));
elemArea = sum(dvolu);
solNorm = Ue(:,1)'*(dvolu.*Ue(:,1)) + Ue(:,2)'*(dvolu.*Ue(:,2));

function [gxy,n,convergence] = recomputeGaussPoint(Xe_M,xyg_g,nDeg,coordRef,gxy0,n)
ngauss = size(n,1);
res_x = 1;
res_fun = 1;
gxy = gxy0;
tol = 1e-10;
max_iter = 10;
convergence = true;
iter = 1;
while max(abs(res_fun(:)))>tol && max(abs(res_x(:)))>tol && iter<max_iter
    ac = n(:,:,2)*Xe_M; % ng x 2
    bd = n(:,:,3)*Xe_M;
    invDetJ = reshape(1./(ac(:,1).*bd(:,2)-bd(:,1).*ac(:,2)),[1 1 ngauss]);
    ac = permute(ac,[2 3 1]);
    bd = permute(bd,[2 3 1]);
    invJ = bsxfun(@times,([bd(2,1,:) -bd(1,1,:); -ac(2,1,:) ac(1,1,:) ]),invDetJ); % 2 x 2 x ng
    aux = reshape(permute(xyg_g-n(:,:,1)*Xe_M,[2 1]),[ 1 2 ngauss]); % 1 x 2 x ng
    gxy = gxy0 + reshape(sum(bsxfun(@times,invJ,aux),2),[2 ngauss]);
    V = Vandermonde_LP(nDeg,coordRef);
    invV = inv(V');
    [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(gxy',nDeg);
    N = (invV*[p,dp_dxi,dp_deta])';
    n(:,:,1) = N(1:ngauss,:);
    n(:,:,2) = N(ngauss+1:2*ngauss,:);
    n(:,:,3) = N(2*ngauss+1:end,:);
    res_fun = xyg_g - n(:,:,1)*Xe_M;
    res_x = gxy - gxy0;
    gxy0 = gxy;
    iter = iter +1;
end
if iter == max_iter
    convergence = false;
end
