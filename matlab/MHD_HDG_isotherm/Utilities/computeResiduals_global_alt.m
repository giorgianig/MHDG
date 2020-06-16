function res_face  = computeResiduals_global_alt...
    (Ne,Nfaces,Nv,nv,Fcon,u_tilde,u,Df,Ef,Lf,Hf,fH,L,L0,U,U0,Fdir)

neq = 2;
res_face = zeros(neq*nv,Nfaces);

for ielem = 1:Ne

    % element faces
    Fe = Fcon(ielem,:);

    % indices
    indu = (ielem-1)*neq*Nv + (1:neq*Nv);

    for iface = 1:size(Fe,2)

        % indices
        indf = (iface-1)*neq*nv + (1:neq*nv);
        ind_utilde = (Fe(iface)-1)*neq*nv + (1:neq*nv);
        uef_tilde = u_tilde(ind_utilde);
        
        
        ind_utilde_el = col(bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)'));

        % face matrices
        Ue = U(:,:,ielem);
        U0e = U0(:,ielem);
        Le = L(:,:,ielem);
        L0e = L0(:,ielem);
        Dff = Df(indf,:,ielem);
        Eff = Ef(indf,indf,ielem);
        Lff = Lf(indf,:,ielem);
        Hff = Hf(indf,indf,ielem);
        fHf = fH(indf,ielem); 

        if ~Fdir(ielem,iface)
            % residuals
            res_face(:,Fe(iface)) = res_face(:,Fe(iface)) + ...
                (  Dff*(Ue*u_tilde(ind_utilde_el) + U0e) + ...
                (Hff - Eff)*uef_tilde) - fHf ...
                -Lff*(Le*u_tilde(ind_utilde_el)+L0e);
        end
    end
end

