%% setNeumannFaces

% set faces to Neumann in all the boundary a part of the upper
tol = 1e-8;
if ~exist('Tb_Diriclet')
    error('Impose Neumann faces only in place of Dirichlet')
end
icount = 1;
nf = size(elementFaceInfo.Diriclet,1);
for iface = 1:nf
    
    el = elementFaceInfo.Diriclet(iface,1);
    nodes = elemInfo.faceNodes(elementFaceInfo.Diriclet(iface,2),:);
    if all( abs( X(T(el,nodes),2)-0) > tol )
        elementFaceInfo.Neumann(icount,:) = elementFaceInfo.Diriclet(iface,:);
        Tb_Neumann(icount,:) = Tb_Diriclet(iface,:);
        elementFaceInfo.Diriclet(iface,:) = -1;
        Tb_Diriclet(iface,:) = -1;
        icount = icount + 1;
    end    
end
Tb_Diriclet = reshape(Tb_Diriclet(Tb_Diriclet~=-1),[nf-icount+1,numel(nodes)]);
elementFaceInfo.Diriclet = reshape(elementFaceInfo.Diriclet(elementFaceInfo.Diriclet~=-1),[nf-icount+1,2]);