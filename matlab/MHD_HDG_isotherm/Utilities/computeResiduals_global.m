function res_face  = computeResiduals_global...
    (Ne,Nfaces,Nv,nv,Fcon,u_tilde,u,Df,Ef,Lf,Hf,fH,Fdir)

global Mesh refEl
neq = 2;
res_face = zeros(neq*nv,Nfaces);

for ielem = 1:Ne

    % element faces
    Fe = Fcon(ielem,:);

    % indices
    indu = (ielem-1)*neq*Nv + (1:neq*Nv);

    % elemental solutions
    ue = u(indu);
    
    Xe = Mesh.X(Mesh.T(ielem,:),:);
    [uex,uxex,uyex] = analyticalSolution(Xe);
    b = defineMagneticField(Xe);
    GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1])]);
    GradPerpU = zeros(size(GradU));
    for ip = 1:size(Xe,1)
        bGradU = GradU(:,:,ip)*transpose(b(ip,:));
        GradPerpU(:,:,ip) = GradU(:,:,ip)-bGradU*b(ip,:);
    end
    
    GradPerpU = col(permute(GradPerpU,[2 1 3]));
    for iface = 1:size(Fe,2)

        % indices
        indf = (iface-1)*neq*nv + (1:neq*nv);
        ind_utilde = (Fe(iface)-1)*neq*nv + (1:neq*nv);       
        
        uef_tilde = u_tilde(ind_utilde);

        % face matrices
        Dff = Df(indf,:,ielem);
        Eff = Ef(indf,indf,ielem);
        Lff = Lf(indf,:,ielem);
        Hff = Hf(indf,indf,ielem);
        fHf = fH(indf,ielem);
        
        if ~Fdir(ielem,iface)
            % residuals
            res_face(:,Fe(iface)) = res_face(:,Fe(iface)) + ...
                (  Dff*ue + (Hff - Eff)*uef_tilde -fHf ...
                -Lff*GradPerpU);
        end
    end
end

