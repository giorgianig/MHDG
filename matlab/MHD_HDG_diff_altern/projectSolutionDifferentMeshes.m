function u0 = projectSolutionDifferentMeshes(meshStart,solStart)

global neq Mesh refEl axisym
X = Mesh.X;
T = Mesh.T;
lscale = 1;%Mesh.lscale;
Sstart  = load(solStart);
MStart = load(meshStart);

% if axisym && min(MStart.X(:,1))<0
%     % apply translation in x to get min(X(:,1)) = 1
%    X(:,1) = X(:,1) - min(X(:,1))+1;
% end

refElstart = createReferenceElement(1,size(MStart.T,2),[]);
u = transpose(reshape(Sstart.u0,neq,numel(MStart.T)));
u_proj = evalDGapproximationAtPoints(X,u,MStart.X/lscale,MStart.T,refElstart);
% disp(' ')
% aux_Gam = evalDGapproximationAtPoints(X,Sstart.u0(2:2:end),MStart.X/lscale,MStart.T,refElstart);

Ne = size(T,1);
Nv = size(T,2);
u0 = zeros(Nv*Ne,neq);
for iel = 1:size(T,1)
    ind = (iel-1)*Nv+(1:Nv);
    Te = T(iel,:);
    u0(ind,:) = u_proj(Te,:);
end