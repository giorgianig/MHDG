function [shapeFunctions,gaussWeights,gaussPoints]=...
    ProvaComputeShapeFunctionsReferenceElementQuads(nDeg,coord1d,nOfGaussPoints)


nOfNodes = (nDeg+1)^2;

[z,w] = gaussLegendre(nOfGaussPoints,-1,1);
nIP = length(w); %number of integration points in each direction
nOfGauss = nIP^2;

%Vandermonde matrix
V = Vandermonde_LP2D_Quad(nDeg,coord1d);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,3);
gaussWeights = zeros(nOfGauss,1);
gaussPoints = zeros(nOfGauss,2);

iGauss = 1;
%Integration over [-1,1]^2
for j = 1:nIP
    y = z(j);
    py = orthopoly1D(y,nDeg);
    for i = 1:nIP
%         x = [z(i),z(j)]; % (r,s) coordinates
%         ind = (j-1)*N+i;
        x = z(i);
        px = orthopoly1D(x,nDeg);
        p = col(px*py');
%         [p,p_xi,p_eta] = orthopoly2D_deriv_rst_qua(x,nDeg);
%         N = U\(L\(P*[p,p_xi,p_eta]));
        N = U\(L\(P*p));
        shapeFunctions(:,iGauss,1) = N(:,1)';
%         shapeFunctions(:,iGauss,2) = N(:,2)';
%         shapeFunctions(:,iGauss,3) = N(:,3)';
%         gaussWeights(iGauss) =(w(i)*w(j))*(1-x(2))/2;
        % only for PFEM
%         r = x(1); s = x(2);
%         xi = (1+r)*(1-s)/2-1;
%         gaussPoints(iGauss,:) = [xi, s];
        iGauss = iGauss + 1;
    end
end