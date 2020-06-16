function gradU = calculateElementByElementGradient(u_tilde,F,L,L0)

global neq Mesh   
nd = size(Mesh.X,2);
Ne = size(F,1);            % number of elements
nf = size(F,2);            % number of faces per element
Nv = size(L,1)/neq/nd;        % nodes per element
nv = size(L,2)/(neq*nf);          % nodes per face

% allocation
gradU = zeros(neq*nd*Nv*Ne,1);

for ielem = 1:Ne
    
    % element faces
    Fe = F(ielem,:);
    
    % indices
    ind_u_tilde = reshape(bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)'),neq*nf*nv,1);
    ind_L = (ielem-1)*neq*Nv*nd + (1:neq*Nv*nd);
    % elemental solutions
    gradU(ind_L) = L(:,:,ielem)*u_tilde(ind_u_tilde) + L0(:,ielem);
end
