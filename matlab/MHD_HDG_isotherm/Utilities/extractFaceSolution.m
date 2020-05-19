function u_tilde = extractFaceSolution(u,F,F_dir,refEl)

neq = 2;
Ne = size(F,1);
maxF = max(F(:));
dirFaces = sum(F_dir(:));
unkF = maxF-dirFaces;
Nv = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
u_tilde = zeros(nv*unkF,neq);
alreadyDone = false(unkF,1);

% resize
u = transpose(reshape(u,neq,Ne*Nv));

for iElem = 1:Ne
    for iface = 1:size(refEl.faceNodes,1)
        if ~F_dir(iElem,iface) && ~alreadyDone(F(iElem,iface))
            faceNodes = refEl.faceNodes(iface,:);
            iFace = F(iElem,iface);
            ind_utilde = (iFace-1)*nv +(1:nv);
            ind_u = (iElem-1)*Nv + (1:Nv);
            u_tilde(ind_utilde,:) = u(ind_u(faceNodes),:);
            alreadyDone(F(iElem,iface)) = true;
        end
    end
end

u_tilde = reshape(transpose(u_tilde),neq*nv*unkF,1);