function u = calculateElementByElementSolution(u_tilde,F,U,U0)

global neq                  % number of equations
Ne = size(F,1);            % number of elements
nf = size(F,2);            % number of faces per element
Nv = size(U,1)/neq;        % nodes per element
nv = size(U,2)/(neq*nf);          % nodes per face

% allocation
u = zeros(neq*Nv*Ne,1);

for ielem = 1:Ne
    
    % element faces
    Fe = F(ielem,:);
    
    % indices
    ind_u_tilde = reshape(bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)'),neq*nf*nv,1);
    ind_u = (ielem-1)*neq*Nv + (1:neq*Nv);
    
    % elemental solutions
    u(ind_u) = U(:,:,ielem)*u_tilde(ind_u_tilde) + U0(:,ielem);
end
