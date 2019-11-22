function res = computeError_3D(u,X,T,theta,refElPol,refElTor,ntor)

u = transpose(reshape(u,4,numel(u)/4));
N2d = size(T,1);                                               % Number of 2d elements
Nel = N2d*ntor;                                               % Number of 3d elements
Np2d      = size(refElPol.NodesCoord,1);              % Number of points for 2d element
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Npel = Np2d*Np1dTor;                                       % Number of points for 3d element
error = zeros(1,4);
tdiv = linspace(0,theta,ntor+1);
volu = 0;
solint2 = error;
%% Loop in theta divisions
for itor = 1:ntor
    tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
    
    %% Loop in 2d elements
    for iel = 1:N2d
        Te = T(iel,:);
        Xe = X(Te,:);
        iElem = (itor-1)*N2d+iel;       % 3d numbering of the element
        ind     = (iElem-1)*Npel+1:iElem*Npel;
        ue = u(ind,:);
        
         [elerr,elvolu,solint2e]= computeElementalError(Xe,tel,refElPol,refElTor,ue);
         error  = error + elerr;
         volu(iElem) = elvolu;
         solint2 = solint2+solint2e;
    end
end
% res = sqrt(sum(error.^2))/sum(volu);
res = sqrt(error.^2./solint2);





%_______________________________________________________________________
function [elemL2Norm,volu,solint2] = computeElementalError(Xe,theta,refElPol,refElTor,ue)

global axisym
%Information of the reference element
ipw2d  = refElPol.IPweights;
ipw1dt = refElTor.IPweights1d;
N2d       = refElPol.N;
N1dTor  = refElTor.N1d;
N2dxi    = refElPol.Nxi;
N2deta  = refElPol.Neta;
N1dxTor= refElTor.N1dxi;

% Delta theta
htheta   = theta(end)-theta(1);

% Number of Gauss points in the interior
ngausspol = length(ipw2d);
ngausstor = length(ipw1dt);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

%Compute elemental L2 Norm
elemL2Norm = zeros(1,4);

volu = 0;
solint2 = elemL2Norm;
% Loop in 1d Gauss points
for igtor = 1:ngausstor
    
    N1_g   = N1dTor(igtor,:);                             % toroidal shape functions
    tg       = N1dTor(igtor,:)*theta;
    dvolu1d = 0.5*ipw1dt(igtor)*htheta;           % toroidal contribution to elemental volume
    
    % Loop in 2D Gauss points
    for igpol = 1:ngausspol
        
        N2_g = N2d(igpol,:);
        N2dxi_g = N2dxi(igpol,:);
        N2deta_g = N2deta(igpol,:);
        
        % 3d shape functions
        Ni  = colt(N2_g'*N1_g);
        
        % Solution at Gauss points
        ueg = Ni*ue;
        
        % Jacobian
        J = [N2dxi_g*xe	    N2dxi_g*ye
            N2deta_g*xe   N2deta_g*ye];
        
        % gauss point position
        xg = N2_g*xe;
        yg = N2_g*ye;
        
        % Analytical solution at Gauss points
        u0g = analyticalSolution_3D([xg,yg],tg);
        
        % Error at Gauss points
        errg = ueg-u0g;
        
        % Integration weight
        dvolu=ipw2d(igpol)*det(J)*dvolu1d;
        if axisym
            dvolu = dvolu*xg;
        end
        
        %Contribution of the current integration point to the elemental L2 Norm
        elemL2Norm = elemL2Norm +errg.^2*dvolu;
        volu = volu+dvolu;
        solint2 = solint2+u0g.^2*dvolu;
    end
end
elemL2Norm = sqrt(elemL2Norm);
% elemL2Norm
