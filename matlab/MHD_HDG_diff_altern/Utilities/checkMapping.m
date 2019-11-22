function checkMapping(LL,LL0,UU,UU0,u_tilde,F,T,X)

global refEl
neq = 2;                   % number of equations
Ne = size(F,1);            % number of elements
nf = size(F,2);            % number of faces per element
nv = size(UU,2)/(neq*nf);          % nodes per face

tol = 1e-9;
figure
for ielem = 1:Ne
    
    % element faces
    Fe = F(ielem,:);
    Te = T(ielem,:);
    Xe = X(Te,:);
    % indices
    ind_u_tilde = reshape(bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)'),...
        neq*nf*nv,1);
    ue_tilde = u_tilde(ind_u_tilde) ;

    [uex,uxex,uyex] = analyticalSolution(Xe);
    GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1])]);
    GradU = col(permute(GradU,[2 1 3]));

    % elemental solutions
    u = UU(:,:,ielem)*ue_tilde + UU0(:,ielem);
    L = LL(:,:,ielem)*ue_tilde + LL0(:,ielem);
    uex = col(transpose(analyticalSolution(Xe)));
%     if max(abs(u-uex))>tol, error('Problema di mapping'), end
%     if max(abs(L-GradU))>tol, error('Problema di mapping'), end
    max(abs(u-uex))
    max(abs(L-GradU))
%     plotSolution(Xe,1:size(T,2),L(4:4:end),refEl,30)
%     hold on
end

stop

