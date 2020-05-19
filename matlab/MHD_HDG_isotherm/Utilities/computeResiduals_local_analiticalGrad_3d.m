function residual = computeResiduals_local_analiticalGrad_3d...
    (Fcon,F_dir,flipFace,u_tilde,u0,A,G,Cv,D,E,Edir,H,Hdir,force,...
    B,C,C_dir,L,P,Pb,Q,Qb,invL,refEl,infoFaces)



global Mesh diff_n diff_u theta
global ntor neq refElTor             % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
Nv     = Np1dTor*Np2d;
N2d    = size(Fcon,1);                    % number of elements
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;
Ne = N2d*ntor;
Nf = max(max(Fcon));
nDirFaces =  sum(sum(F_dir));


nd = 3;
X = Mesh.X;
T = Mesh.T;
residual = zeros(2,Ne);


% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
indu_tilde = zeros(Nfp*neq,1);
tdiv = linspace(0,theta,ntor+1);



for itor = 1:ntor
    tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
    for iel = 1:N2d
        
        % element faces
        Fcone = Fcon(iel,:);
        
        iElem = (itor-1)*N2d + iel;
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*[(iel-1)*Np2d,N2d*Np2d+...
            (Fcone-1)*Nfl,N2d*Np2d+(Nf-nDirFaces)*Nfl+(iel-1)*Np2d];
        if itor==ntor
            delta(end) = 1+(iel-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            indu_tilde(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
        end
        
        
        % indices
%         indL = (iElem-1)*neq*Nv*nd + (1:neq*Nv*nd);
        indu = (iElem-1)*neq*Nv + (1:neq*Nv);
        
        % elemental matrices
        Ae      = A(:,:,iElem);
        Cve     = Cv(:,:,iElem);
        Ge      = G(:,:,iElem);
        De      = D(:,:,iElem);
        Ee      = E(:,:,iElem);
        Edire   = Edir(:,iElem);
        He      = H(:,:,iElem);
        Hdire   = Hdir(:,iElem);
        f       = force(:,iElem);
        Be      = B(:,:,iElem);
        Ce      = C(:,:,iElem);
        C_dire  = C_dir(:,iElem);
        Le      = L(:,:,iElem);
        Pe      = P(:,:,iElem);
        Qe      = Q(:,:,iElem);
        Pbe     = Pb(:,:,iElem);
        Qbe     = Qb(:,:,iElem);
        invLe   = invL(:,:,iElem);
        
        % elemental solutions
        ue = u0(indu);
        ue_tilde = u_tilde(indu_tilde) ;
        
        Xe = X(T(iel,:),:);
        [uex,uxex,uyex,utex] = analyticalSolution3d(Xe,tel);
        uex = col(transpose(uex));

        
        GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1]),permute(utex,[2 3 1])]);
        GradUe = col(permute(GradU,[2 1 3]));
        ue_tilde_ex = zeros(Nfp*neq,1);
        for iface = 1:nf+2
            if iface==1
                ind_ff = 1:Np2d*neq;                                 % assembly face-face
                ind_fe = 1:Np2d*neq;                                 % assembly face-element (for var)
            elseif iface==nf+2
                ind_ff   = Np2d*neq+nf*Nfl*neq + (1:Np2d*neq);
                ind_fe   = Np2d*(Np1dTor-1)*neq + (1:Np2d*neq);
            else
                nodesv      = refElTor.faceNodes3(iface-1,:);
                ind_fe      = col(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
                ind_ff      = Np2d*neq + (iface-2)*Np1dPol*Np1dTor*neq+ (1:Np1dPol*Np1dTor*neq);
            end
            ue_tilde_ex(ind_ff) = uex(ind_fe);
        end
        
        
        flipFace_e = flipFace(iel,:);
        for iface = 1:nf
            if flipFace_e(iface)
                ind_v_L = Np2d*neq + (iface-1)*Nfl*neq+ (1:Nfl*neq);
                ue_tilde(ind_v_L) = ue_tilde(ind_v_L(perm));
            end
        end
        diffdiag = repmat(col(repmat([diff_n diff_u],nd,1)),Nv,1);
        Diff = diag(diffdiag);
        
        
        % residuals
        residual(1,iElem) = norm( Le*GradUe+Be*ue-Ce*ue_tilde-C_dire );
        residual(2,iElem) = norm( (-Cve+De+Ge-((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*Be)*ue+...
            (He-Ee+((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*Ce)*ue_tilde+...
            Hdire-Edire -f+((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*C_dire ...
            );
        %      stop
    end
end


