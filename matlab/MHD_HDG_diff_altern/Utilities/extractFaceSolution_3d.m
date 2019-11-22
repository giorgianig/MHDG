function u_tilde = extractFaceSolution_3d(u,F,F_dir,infoFaces,refEl)

global ntor refElTor

neq  = 2;
N2D  = size(F,1);
Np2D = size(refEl.NodesCoord,1);
Ne   = N2D*ntor;
Nf   = max(max(F));
Np1Dpol = size(refEl.NodesCoord1d,1);
Np1Dtor = size(refElTor.NodesCoord1d,1);
Nfl  = Np1Dpol*Np1Dtor;
Nfdir = sum(sum(F_dir));
unkF = Nf-Nfdir;
Np   = Np2D*Np1Dtor;
nut  = ntor*(Nfl*unkF + Np2D*N2D);

indl = 1:Nfl;
ind2 = 1:Np2D;
ind3 = 1:Np;
% resize
u_tilde = zeros(nut,neq);
u = transpose(reshape(u,neq,Ne*Np));
c = 0;

% exterior faces infos
bnames = fieldnames(infoFaces);

for itor = 1:ntor
    % Poloidal faces
    for iel = 1:N2D
        iElem = (itor-1)*N2D + iel;
        ind_ue = (iElem-1)*Np  + ind3;
        u_tilde(c+ind2,:) = u(ind_ue(ind2),:);
        c = c+Np2D;
    end
    
    % Toroidal interior faces
    for iFace = 1:size(infoFaces.interiorFaces,1)
        iel = infoFaces.interiorFaces(iFace,1);
        ifa = infoFaces.interiorFaces(iFace,2);
        iElem = (itor-1)*N2D + iel;
        ind_ue = (iElem-1)*Np  + ind3;
        u_tilde(c+indl,:) = u(ind_ue(refElTor.faceNodes3(ifa,:)),:);
        c = c+Nfl;
    end
    
    % Toroidal exterior faces
    for ib = 1:numel(bnames)
        nm = bnames{ib};
        if ~strcmpi(nm(1:3),'ext'),continue,end
        name = nm(15:end);
        faces = infoFaces.(['exteriorFaces_', name]);
        for iFace = 1:size(faces,1)
            iel = faces(iFace,1);
            ifa = faces(iFace,2);
            if ( ~F_dir(iel,ifa) )
                iElem = (itor-1)*N2D + iel;
                ind_ue = (iElem-1)*Np  + ind3;
                u_tilde(c+indl,:) = u(ind_ue(refElTor.faceNodes3(ifa,:)),:);
                c = c+Nfl;
            end
        end
    end
end

u_tilde = reshape(transpose(u_tilde),nut*neq,1);