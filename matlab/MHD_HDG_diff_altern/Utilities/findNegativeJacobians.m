% ONLY WORKS WITH NON-CURVED ELEMENTS
clear all
home
close all
% addpath('../Saves/BarcelonaCoarseData/h10/')
load ../../Meshes/KT_profile/KT_profile_P5.mat

nOfElements = size(T,1);
negative = false(size(T,1),1);
small_jac = false(size(T,1),1);
% loop in elements
for ielem =  1:nOfElements
    Te = T(ielem,1:3);
    Xe = X(Te,:);
    v1 = Xe(1,:);  v2 = Xe(2,:);  v3 = Xe(3,:);
    J = [(v2-v1)/2 ; (v3-v1)/2];
    detJ = det(J);
    if detJ<1e-3
        if detJ<0
        negative(ielem) = true;
        end
        small_jac(ielem) = true;
    end
end

if any(negative)
    plotMesh(X,T);
    hold on
    patch('Faces',T(negative,:),'Vertices',X,'FaceVertexCData',1,...
        'FaceColor','flat','EdgeAlpha',0);
else
    disp('No negative jacobians')
    if any(small_jac)
        disp([num2str(nnz(small_jac)) ' elements with detJ<1e-3'])
    else
        disp('No small jacobians')
    end 
end



