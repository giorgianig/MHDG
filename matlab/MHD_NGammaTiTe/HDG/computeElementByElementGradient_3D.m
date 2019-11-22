function q = computeElementByElementGradient_3D(ntor,lambda,QQ,Q0,T,F,Fdir,refElPol,refElTor)

global neq
Np2d      = size(refElPol.NodesCoord,1);              % Number of points for 2d element
Np1dPol   = size(refElPol.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor   = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Npel      = Np2d*Np1dTor;                                       % Number of points for 3d element
N2d       = size(T,1);                                                % Number of 2d elements
Ne        = N2d*ntor;                                               % Number of 3d elements
nf        = size(refElPol.faceNodes,1);                      % Number of faces in the 2d element
Ndim      = 3;                                                            % Number of spatial dimensions
Nfl       = Np1dPol*Np1dTor;                                    % Number of nodes for lateral faces
Nfp       = Np2d*2+nf*Np1dPol*Np1dTor;                  % Number of nodes in all the 3d faces of 1 element
Nf        = max(max(F));                                              % Number of faces in the 2d plane
nDirFaces =  sum(sum(Fdir));


q = zeros(Ne*Ndim*Npel*neq,1);
ind_fac = zeros(Nfp*neq,1);
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
for itor = 1:ntor
    for iel = 1:N2d
        iElem = (itor-1)*N2d+iel; % 3d numbering of the element
        ind_gra = (iElem-1)*Npel*Ndim*neq + (1:Npel*Ndim*neq);
        Fe = F(iel,:);
        
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+...
            neq*[(iel-1)*Np2d,N2d*Np2d+(Fe-1)*Nfl,N2d*Np2d+...
            (Nf-nDirFaces)*Nfl+(iel-1)*Np2d];
        
        if itor==ntor
            delta(end) = 1+(iel-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            ind_fac(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
        end
        
         q(ind_gra) = QQ(:,:,iElem)*lambda(ind_fac)+Q0(:,iElem);
    end
end

