function L_gradv = calculateElementByElementVelocityGradient(u_tilde,ro,F,L,L0,Lro,Lf,U)

Ne = size(F,1);
Nv = size(U,1)/2;          % velocity nodes per element
nv = size(U,2)/6;          % velocity nodes per face

% allocation
L_gradv = zeros(4*Nv*Ne,1);

for ielem = 1:Ne
    
    % element faces
    Fe = F(ielem,:);
    
    % indices
    ind_u_tilde = reshape(bsxfun(@plus,(Fe-1)*2*nv,(1:2*nv)'),6*nv,1);
    ind_L = (ielem-1)*4*Nv + (1:4*Nv);
    
    % elemental solutions
    L_gradv(ind_L) = L(:,:,ielem)*u_tilde(ind_u_tilde) + Lro(:,:,ielem)*ro(ielem) + ...
        L0(:,:,ielem) + Lf(:,:,ielem);

end