function checkMapping3d(LL,LL0,UU,UU0,u_tilde,F,F_dir,T,X)

global refEl theta ntor
global neq refElTor             % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
% Nv     = Np1dTor*Np2d;
N2d    = size(F,1);                    % number of elements
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;
Ne = N2d*ntor;
Nf = max(max(F));
nDirFaces =  sum(sum(F_dir));
nf = size(F,2);            % number of faces per element

% tol = 1e-9;
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
indu_tilde = zeros(Nfp*neq,1);
tdiv = linspace(0,theta,ntor+1);

for itor = 1:ntor
    tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
    for ielem = 1:N2d
        
        % element faces
        Fe = F(ielem,:);
        Te = T(ielem,:);
        Xe = X(Te,:);

        iElem = (itor-1)*N2d + ielem;
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*[(ielem-1)*Np2d,N2d*Np2d+...
            (Fe-1)*Nfl,N2d*Np2d+(Nf-nDirFaces)*Nfl+(ielem-1)*Np2d];
        if itor==ntor
            delta(end) = 1+(ielem-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            indu_tilde(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
        end
        
        
        % indices
%         indL = (iElem-1)*neq*Nv*nd + (1:neq*Nv*nd);
%         indu = (iElem-1)*neq*Nv + (1:neq*Nv);

        ue_tilde = u_tilde(indu_tilde) ;
        
        [uex,uxex,uyex,utex] = analyticalSolution3d(Xe,tel);
        GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1]),permute(utex,[2 3 1])]);
        GradU = col(permute(GradU,[2 1 3]));
        
        % elemental solutions
        u = UU(:,:,iElem)*ue_tilde + UU0(:,iElem);
        L = LL(:,:,iElem)*ue_tilde + LL0(:,iElem);
        uex = col(transpose(uex));
        %     if max(abs(u-uex))>tol, error('Problema di mapping'), end
        %     if max(abs(L-GradU))>tol, error('Problema di mapping'), end
        max(abs(u-uex))
        max(abs(L-GradU))
        %     plotSolution(Xe,1:size(T,2),L(4:4:end),refEl,30)
        %     hold on
    end
end
stop

