function plotFaceVelocity(X,T,F,u_tilde,refElv)

nElems = size(T,1);
nFaces = max(max(F));
faceNodes = refElv.faceNodes;
nv = size(faceNodes,2);
check = false(nFaces,1);

% initialize figure
figure,clf
plotMesh(X,T)
box on; axis equal tight;
hold on

% loop in elements
for ielem = 1:nElems

    Fe = F(ielem,:);
    Te = T(ielem,:);

    for iface = 1:3

        face = Fe(iface);

        if ~check(face)

            Tf = Te(faceNodes(iface,:));
            xf = X(Tf,1);
            yf = X(Tf,2);
            
            % index
            indx = (face-1)*2*nv + (1:2:2*nv-1);
            indy = (face-1)*2*nv + (2:2:2*nv);

            % plot
            quiver(xf,yf,u_tilde(indx),u_tilde(indy),0.5);
            hold on
            
            % actualize check
            check(face) = true;
        end
    end
end
axis equal; axis([-0.5,1.5,-0.5,1.5]);