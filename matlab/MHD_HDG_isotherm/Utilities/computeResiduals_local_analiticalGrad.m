function residual = computeResiduals_local_analiticalGrad...
    (Ne,Nv,nv,Fcon,flipFace,u_tilde,u0,A,G,Cv,D,E,Edir,H,Hdir,force,...
    B,C,C_dir,L,P,Pb,Q,Qb,dt,invL)
global Mesh diff_n diff_u

X = Mesh.X;
T = Mesh.T;


nf  = 3;  % nOfFaces for HexaRefElem !!!
neq = 2;
nd = 2;
residual = zeros(2,Ne);
% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

for iface = 1:nf
    ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

for ielem = 1:Ne
    
    % element faces
    Fcone = Fcon(ielem,:);
    
    % indices
    indL = (ielem-1)*nd*neq*Nv + (1:nd*neq*Nv);
    indu = (ielem-1)*neq*Nv + (1:neq*Nv);
    indu_tilde = reshape(bsxfun(@plus,(Fcone-1)*neq*nv,(1:neq*nv)'),neq*nf*nv,1);
    
    % elemental matrices 
    Ae      = A(:,:,ielem);
    Cve      = Cv(:,:,ielem);
    Ge      = G(:,:,ielem);
    De      = D(:,:,ielem);
    Ee      = E(:,:,ielem);
    Edire   = Edir(:,ielem);
    He      = H(:,:,ielem);
    Hdire   = Hdir(:,ielem);
    f       = force(:,ielem);
    Be      = B(:,:,ielem);
    Ce      = C(:,:,ielem);
    C_dire  = C_dir(:,ielem);
    Le      = L(:,:,ielem);
    Pe      = P(:,:,ielem);
    Qe      = Q(:,:,ielem);
    Pbe     = Pb(:,:,ielem);
    Qbe     = Qb(:,:,ielem);
    invLe   = invL(:,:,ielem);
    
    % elemental solutions
    ue = u0(indu);
    ue_tilde = u_tilde(indu_tilde) ;
    
    Xe = X(T(ielem,:),:);
    [uex,uxex,uyex] = analyticalSolution(Xe);
    GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1])]);
    GradUe = col(permute(GradU,[2 1 3]));
     
    uex = col(transpose(uex));
    
    flipFace_e = flipFace(ielem,:);
    for iface = 1:nf
        if flipFace_e(iface)
            ue_tilde(ind_v_L(iface,:)) = ue_tilde(ind_v_L(iface,perm));
        end
    end
    diffdiag = repmat(col(repmat([diff_n diff_u],nd,1)),Nv,1);
    Diff = diag(diffdiag);
    
    % residuals
    residual(1,ielem) = norm( Le*GradUe+Be*ue-Ce*ue_tilde-C_dire );
    residual(2,ielem) = norm( (-Cve+De+Ge-((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*Be)*ue+...
                              (He-Ee+((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*Ce)*ue_tilde+...
                              Hdire-Edire -f+((Pe-Pbe)-(Qe-Qbe))*Diff*invLe*C_dire ...
                              );   
%      stop
end
    


