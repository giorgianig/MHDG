function upol = extractSolutionInAPoloidalPlane(u,T,refElPol,refElTor,iplane)

global ntor

if iplane<1 && iplane>ntor
    error('Plane not valid')
end

Np2d      = size(refElPol.NodesCoord,1);              % Number of points for 2d element
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nv     = Np1dTor*Np2d;                                      % Number of points for 3d element
N2d   = size(T,1);                                                % Number of 2d elements

upol = zeros(N2d*Np2d,1);
for iel = 1:N2d
        iElem = (iplane-1)*N2d+iel; % 3d numbering of the element
        ind_ue = (iElem-1)*Nv + (1:Nv);
        ind_upol = (iel-1)*Np2d + (1:Np2d);
        ue = u(ind_ue);
        upol(ind_upol) = ue(1:Np2d);
end

