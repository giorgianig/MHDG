function [error eg domArea] = calculateL2ErrorTwoSolutionNested(X1,T1,u1,refEl1,X2,T2,u2,refEl2)

% Compute the L2 error element by element for two nested meshes. The
% routine automatically detect which one is the reference solution (that is
% the one with more elements). Works until 4 nesting levels.
% In the following:  M is the Mother mesh, C is the Child mesh (nested)
% 

if size(T1,1) > size(T2,1)
    T_C = T1;
    X_C = X1;
    u_C = u1;
    refEl_C = refEl1;
    T_M = T2;
    X_M = X2;
    u_M = u2;
    refEl_M = refEl2;
else
    T_C = T2;
    X_C = X2;
    u_C = u2;
    refEl_C = refEl2;
    T_M = T1;
    X_M = X1;
    u_M = u1;
    refEl_M = refEl1;
end

% Mesh info
nElem_C = size(T_C,1);
nElem_M = size(T_M,1);
nestLev = floor(nElem_C/nElem_M);

% elements to compute
elements = 1:nElem_M;
enlacMap_C = 1:nElem_C;

% cross connectivity
[C Rot] = createCrossConnectivity(X_M,T_M,X_C,T_C,refEl_M.faceNodes);

% number of nodes per element
nOfNodes_M = size(refEl_M.NodesCoord,1); % N of nodes mother mesh
nOfNodes_C = size(refEl_C.NodesCoord,1); % N of nodes child mesh

% reshape velocity
u_M = transpose(reshape(u_M,2,nElem_M*nOfNodes_M));
u_C = transpose(reshape(u_C,2,nElem_C*nOfNodes_C));

% Compute shape functions at Gauss points of the child mesh
coordGauss = refEl_C.IPcoordinates;  % Gauss points child mesh
coordRef = refEl_M.NodesCoord;       % Nodes on the mother mesh
nOfGauss = size(coordGauss,1);
nDeg = refEl_M.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');
shapeFunctions_cell = cell(nestLev,1);
coordGauss_map_cell = shapeFunctions_cell;
for iLev = 1:nestLev
    coordGauss_map = TraslateScaleRotate(coordGauss,nestLev,iLev);
    shapeFunctions = zeros(nOfGauss,nOfNodes_M);
    for ipoint = 1:nOfGauss
        [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(coordGauss_map(ipoint,:),nDeg);
        N = (invV*[p,dp_dxi,dp_deta])';
        shapeFunctions(ipoint,:,1) = N(1,:);
        shapeFunctions(ipoint,:,2) = N(2,:);
        shapeFunctions(ipoint,:,3) = N(3,:);
    end
    shapeFunctions_cell{iLev} = shapeFunctions;
    coordGauss_map_cell{iLev} = coordGauss_map;
end

%free memory
nDeg_C = refEl_C.degree;
loop_ind = 0;

% load rotation structure
load rotationStructure.mat

% allocate vectors
error2 = zeros(nElem_M,1);
solNorm = zeros(nElem_M,1);
domArea = ones(nElem_M,1);

% Loop in elements
for ielem = 1:numel(elements)
    loop_ind = loop_ind+1;
    iElem = elements(ielem);
    Ce_M = C(ielem,:);
    Rote_M = Rot(ielem,:);
    Te_M = T_M(iElem,:);
    Xe_M = X_M(Te_M,:);
    ind_M =  (iElem-1)*nOfNodes_M+1:iElem*nOfNodes_M;
    ue_M = u_M(ind_M,:);
    elemError = 0;
    elemArea = 0;
    elemSolNorm = 0;
    for n = 1:nestLev
        jelem = Ce_M(n);
        Te_C = T_C(enlacMap_C(jelem),:);
        Xe_C = X_C(Te_C,:);
        rot = rotationStructure{Rote_M(n)}{nDeg_C};
        ind_C =  (enlacMap_C(jelem)-1)*nOfNodes_C+1:enlacMap_C(jelem)*nOfNodes_C;
        ue_C = u_C(ind_C,:);
        [error_n solNorm_n area_n] = calculateElementalError(ue_C,ue_M,Xe_C,Xe_M,refEl_C,...
            shapeFunctions_cell{n},rot,coordGauss_map_cell{n},nDeg,...
            coordRef);
        elemError =  elemError + error_n;
        elemArea = elemArea + area_n;
        elemSolNorm = elemSolNorm + solNorm_n;
    end
    error2(iElem) = elemError;
    solNorm(iElem) = elemSolNorm;
    domArea(iElem) = elemArea;
end

error = sqrt(error2/sum(solNorm));
eg = sqrt(sum(error2)/sum(solNorm));

function [elemError solNorm elemArea] = calculateElementalError...
    (ue_C,ue_M,Xe_C,Xe_M,referenceElement,n,rot,gauss_map_M,nDeg,coordRef)
% ue_C: reference solution (on the nested mesh)
% ue_M: solution (on the mother mesh)

%Information of the reference element
IPw = referenceElement.IPweights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;

xyg_xi = Nxi*Xe_C(rot,:);
xyg_eta = Neta*Xe_C(rot,:);
detJ = xyg_xi(:,1).*xyg_eta(:,2) - xyg_xi(:,2).*xyg_eta(:,1);
dvolu = detJ.*IPw;
xyg_C = N*Xe_C(rot,:);
xyg_M = n(:,:,1)*Xe_M;
xyg_err = (xyg_C - xyg_M);
if max(abs(xyg_err(:)))>1e-10
    gxy0 = transpose(gauss_map_M);
    [gxy,n,convergence] = recomputeGaussPoint(Xe_M,xyg_C,nDeg,coordRef,gxy0,n);
    if ~convergence
        disp('not converging')
    end
end
Ue = N*ue_C(rot,:);
ue = n(:,:,1)*ue_M;
err_gauss = Ue-ue;
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

function xy_ts =  TraslateScaleRotate(xy,nestLev,nLoc)
% translate and scale the gauss point coordinates
% xy: coordinates of the gauss points in the reference element
% nestLev: level of nested mesh
% nLoc: local number of the nested element

d = 2/(sqrt(nestLev)); % length of the enlaced element edge
if ScaleRotate(nestLev,nLoc) % scale and translate
    A = [-1 0; 0 -1];
    xy_ts = transpose(A*((xy+1)/sqrt(nestLev))' -1 + d);
else           % scale and rotate traslate
    xy_ts = (xy+1)/sqrt(nestLev) -1;
end
% now translate in the correct location given by the local number
[nx, ny] = locateRefEl(nestLev,nLoc);
xy_ts(:,1) = xy_ts(:,1) + nx*d;
xy_ts(:,2) = xy_ts(:,2) + ny*d;

function res = ScaleRotate(nestLev,nLoc)
% 0: scale
% 1: scale & rotate
switch nestLev
    case 4
        switch nLoc
            case {1,3,4}
                res = 0;
            case 2
                res = 1;
        end
    case 16
        switch nLoc
            case {1,3,5,7,8,10,12,13,15,16}
                res = 0;
            case {2,4,6,9,11,14}
                res = 1;
        end
    case 64
        switch nLoc
            case {1,3,5,7,9,11,13,15,16,18,20,22,24,26,28,29,31,33,35,37,39,40,42,...
                    44,46,48,49,51,53,55,56,58,60,61,63,64}
                res = 0;
            case {2,4,6,8,10,12,14,17,19,21,23,25,27,30,32,34,36,36,38,41,43,45,47,50,52,54,57,59,62}
                res = 1;
        end
    case 256
        switch nLoc
            case {1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,...
                    61,63,65,67,69,71,73,75,77,79,81,83,85,87,88,90,92,94,96,98,100,102,104,106,108,110,112,...
                    113,115,117,119,121,123,125,127,129,131,133,135,136,138,140,142,144,146,148,150,152,154,156,...
                    157,159,161,163,165,167,169,171,173,175,176,178,180,182,184,186,188,190,192,193,195,197,199,...
                    201,203,205,207,208,210,212,214,216,218,220,221,223,225,227,229,231,232,234,236,238,240,241,243,...
                    245,247,248,250,252,253,255,256}
                res = 0;
            case {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,33,35,37,39,41,43,45,47,49,51,53,55,57,59,...
                    62,64,66,68,70,72,74,76,78,80,82,84,86,89,91,93,95,97,99,101,103,105,107,109,111,...
                    114,116,118,120,122,124,126,128,130,132,134,137,139,141,143,145,147,149,151,153,155,...
                    158,160,162,164,166,168,170,172,174,177,179,181,183,185,187,189,191,194,196,198,200,202,204,206,...
                    209,211,213,215,217,219,222,224,226,228,230,233,235,237,239,242,244,246,249,251,254}
                res = 1;
        end
end



function [nx, ny] = locateRefEl(nestLev,nLoc)
% map the enlaced reference element in the original one


switch nestLev
    case 4
        % x
        switch nLoc
            case {1,2,4}
                nx = 0;
            case 3
                nx = 1;
        end
        % y
        switch nLoc
            case {1,2,3}
                ny = 0;
            case 4
                ny = 1;
        end
    case 16
        % x
        switch nLoc
            case {1,2,8,9,13,14,16}
                nx = 0;
            case {3,4,10,11,15}
                nx = 1;
            case {5,6,12}
                nx = 2;
            case 7
                nx = 3;
        end
        % y
        switch nLoc
            case {1,2,3,4,5,6,7}
                ny = 0;
            case {8,9,10,11,12}
                ny = 1;
            case {13,14,15}
                ny = 2;
            case 16
                ny = 3;
        end
    case 64
        % x
        switch nLoc
            case {1,2,16,17,29,30,40,41,49,50,56,57,61,62,64}
                nx = 0;
            case {3,4,18,19,31,32,42,43,51,52,58,59,63}
                nx = 1;
            case {5,6,20,21,33,34,44,45,53,54,60}
                nx = 2;
            case {7,8,22,23,35,36,46,47,55}
                nx = 3;
            case {9,10,24,25,37,38,48}
                nx = 4;
            case {11,12,26,27,39}
                nx = 5;
            case {13,14,28}
                nx = 6;
            case 15
                nx = 7;
        end
        % y
        switch nLoc
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
                ny = 0;
            case {16,17,18,19,20,21,22,23,24,25,26,27,28}
                ny = 1;
            case {29,30,31,32,33,34,35,36,37,38,39}
                ny = 2;
            case {40,41,42,43,44,45,46,47,48}
                ny = 3;
            case {49,50,51,52,53,54,55}
                ny = 4;
            case {56,57,58,59,60}
                ny = 5;
            case {61,62,63}
                ny = 6;
            case 64
                ny = 7;
        end
        
    case 256
        % x
        switch nLoc
            case {1,2,32,33,61,62,88,89,113,114,136,137,157,158,176,177,193,194,208,209,221,222,232,233,241,242,...
                    248,249,253,254,256}
                nx = 0;
            case {3,4,34,35,63,64,90,91,115,116,138,139,159,160,178,179,195,196,210,211,223,224,234,235,243,244,...
                    250,251,255}
                nx = 1;
            case {5,6,36,37,65,66,92,93,117,118,140,141,161,162,180,181,197,198,212,213,225,226,236,237,245,246,...
                    252}
                nx = 2;
            case {7,8,38,39,67,68,94,95,119,120,142,143,163,164,182,183,199,200,214,215,227,228,238,239,247}
                nx = 3;
            case {9,10,40,41,69,70,96,97,121,122,144,145,165,166,184,185,201,202,216,217,229,230,240}
                nx = 4;
            case {11,12,42,43,71,72,98,99,123,124,146,147,167,168,186,187,203,204,218,219,231}
                nx = 5;
            case {13,14,44,45,73,74,100,101,125,126,148,149,169,170,188,189,205,206,220}
                nx = 6;
            case {15,16,46,47,75,76,102,103,127,128,150,151,171,172,190,191,207}
                nx = 7;
            case {17,18,48,49,77,78,104,105,129,130,152,153,173,174,192}
                nx = 8;
            case {19,20,50,51,79,80,106,107,131,132,154,155,175}
                nx = 9;
            case {21,22,52,53,81,82,108,109,133,134,156}
                nx = 10;
            case {23,24,54,55,83,84,110,111,135}
                nx = 11;
            case {25,26,56,57,85,86,112}
                nx = 12;
            case {27,28,58,59,87}
                nx = 13;
            case {29,30,60}
                nx = 14;
            case {31}
                nx = 15;
        end
        % y
        switch nLoc
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}
                ny = 0;
            case {32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60}
                ny = 1;
            case {61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87}
                ny = 2;
            case {88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112}
                ny = 3;
            case {113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135}
                ny = 4;
            case {136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156}
                ny = 5;
            case {157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175}
                ny = 6;
            case {176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192}
                ny = 7;
            case {193,194,195,196,197,198,199,200,201,202,203,204,205,206,207}
                ny = 8;
            case {208,209,210,211,212,213,214,215,216,217,218,219,220}
                ny = 9;
            case {221,222,223,224,225,226,227,228,229,230,231}
                ny = 10;
            case {232,233,234,235,236,237,238,239,240}
                ny = 11;
            case {241,242,243,244,245,246,247}
                ny = 12;
            case {248,249,250,251,252}
                ny = 13;
            case {253,254,255}
                ny = 14;
            case 256
                ny = 15;
        end
end
