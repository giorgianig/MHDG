function [C Rot] = createCrossConnectivity(X_1,T_1,X_2,T_2,faceNodes)

% create connectivity of the mother mesh on the enlaced one to compute
% error with a reference mesh
% faceNodes is optional: if used, the whole boundary element is used to
% search for nested elements (to be used for nesting levels > 2)

% find which is the mother and the child
nElem_1 = size(T_1,1);
nElem_2 = size(T_2,1);
if nElem_1 > nElem_2
    %     if mod(nElem_1,nElem_2)
    %         error('The two meshes are not enlaced')
    %     end
    X_C = X_1;
    T_C = T_1;
    X_M = X_2;
    T_M = T_2;
    nElem_C = nElem_1;
    nElem_M = nElem_2;
elseif nElem_1 < nElem_2
    %     if mod(nElem_2,nElem_1)
    %         error('The two meshes are not enlaced')
    %     end
    X_C = X_2;
    T_C = T_2;
    X_M = X_1;
    T_M = T_1;
    nElem_C = nElem_2;
    nElem_M = nElem_1;
else
    error('The two meshes have the same number of elements')
end
clear mesh_1 mesh_2 nElem_1 nElem_2

% define the level of nesting
ref_nestLev = [4 16 64 256];
nestLev = floor(nElem_C/nElem_M);
if ~any(isequal(nestLev,ref_nestLev))
       [v_aux p_aux] = min(abs(nestLev - ref_nestLev));
       nestLev = ref_nestLev(p_aux);
end 

elements_C = 1:size(T_C,1);
C = zeros(nElem_M,nestLev);
Rot = C;
x_C = X_C(:,1);
y_C = X_C(:,2);
aux_T = T_C(elements_C,1:3);
x_vert_C = x_C(aux_T);
y_vert_C = y_C(aux_T);
x_bar_C = sum(x_vert_C,2)/3;
y_bar_C = sum(y_vert_C,2)/3;
for ielem_M = 1:nElem_M
    if isempty(faceNodes)
        coord = X_M(T_M(ielem_M,1:3),:); % use just the vertices
    else
        elemBound = [1 faceNodes(1,:) 2  faceNodes(2,:) 3 faceNodes(3,:)];
        coord = X_M(T_M(ielem_M,elemBound),:); % use the whole elem boundary
    end
    [inElem onElem] = inpolygon(x_bar_C,y_bar_C,coord(:,1),coord(:,2));
    insideElemInd = any([inElem onElem],2);
    insideElem = elements_C(insideElemInd);
    if numel(insideElem)~=nestLev
        error('Something wrong in the element search')
    end
    [insideElemOrd insideElemRot] = findLocalOrder(nestLev,X_M(T_M(ielem_M,1),:),...
        x_vert_C(insideElemInd,:),y_vert_C(insideElemInd,:),T_C(insideElem,1:3));
    C(ielem_M,:) = insideElem(insideElemOrd);
    Rot(ielem_M,:) = insideElemRot;
%     elements_C = elements_C(~insideElemInd);
end

function [insideElemOrd insideElemRot] = findLocalOrder(nestLev,xy,x_C,y_C,T_C)

node = T_C(comp(x_C,xy(1)) & comp(y_C,xy(2)));
insideElemOrd = zeros(nestLev,1);
insideElemRot = insideElemOrd;
row_index = sqrt(nestLev)-1;
ind = 1;
for i = 1:sqrt(nestLev)
    [elem faceNodes node rot] = findFirstRowElem(node,T_C,insideElemOrd(1:ind-1));
    insideElemOrd(ind) = elem;
    insideElemRot(ind) = rot;
    ind = ind+1;
    for j = 1:row_index
        % next element across the hypotenuse
        [elem faceNodes rot] = findSuccessiveElemHYP(elem,faceNodes,T_C);
        insideElemOrd(ind) = elem;
        insideElemRot(ind) = rot;
        ind = ind+1;
        % next element across the cathetus
        [elem faceNodes rot] = findSuccessiveElemCAT(elem,faceNodes,T_C);
        insideElemOrd(ind) = elem;
        insideElemRot(ind) = rot;
        ind = ind+1;
    end
    row_index = row_index-1;
end
% check
if ~isempty(setdiff((1:nestLev),unique(insideElemOrd)))
    error('Something is wrong in the reorder process')
end

function [elem face node rot] = findFirstRowElem(node,T_C,elems)
% find the first element of the row, the face across which continue the
% process and the node to find the successive row

% element
aux_elem = find(any(T_C == node,2));
elem = setdiff(aux_elem,elems);

% face
Te = T_C(elem,:);
face = Te(Te~=node);
% node
aux_node = Te==node;
if all(aux_node == [1 0 0])
    node = Te(3);
    rot = 1;
elseif all(aux_node == [0 1 0])
    node = Te(1);
    rot = 2;
elseif all(aux_node == [0 0 1])
    node = Te(2);
    rot = 3;
end

function [elem face rot] = findSuccessiveElemHYP(felem,faceNodes,T)
% find the successive element across the hypotenuse of the preceeding
% element: felem=preceding element, face=node of the hypotenuse,
% T=connectivity of the children mesh in the element of the mother mesh

index_matrix = (T == faceNodes(1) | T == faceNodes(2));
aux_elem = find(sum(index_matrix,2)==2);
elem = aux_elem(aux_elem~=felem);
Te = T(elem,:);
aux_face = index_matrix(elem,:);
if  all(aux_face == [1 1 0])
    face = Te([2 3]); % face 2
    rot = 3;
elseif all(aux_face ==[0 1 1])
    face = Te([1 3]); % face 3
    rot = 1;
elseif all(aux_face == [1 0 1])
    face = Te([1 2]); % face 1
    rot = 2;
end

function [elem face rot] = findSuccessiveElemCAT(felem,faceNodes,T)
% find the successive element across the cathetus of the preceeding
% element: felem=preceding element, faceNodes= nodes of the cathetus,
% T=connectivity of the children mesh in the element of the mother mesh

index_matrix = (T == faceNodes(1) | T == faceNodes(2));
aux_elem = find(sum(index_matrix,2)==2);
elem = aux_elem(aux_elem~=felem);
Te = T(elem,:);
aux_face = index_matrix(elem,:);
if  all(aux_face == [1 1 0])
    face = Te([1 3]); % face 3
    rot = 2;
elseif all(aux_face ==[0 1 1])
    face = Te([1 2]); % face 1
    rot = 3;
elseif all(aux_face == [1 0 1])
    face = Te([2 3]);% face 2
    rot = 1;
end

function res = comp(x1,x2)
% compare two coordinates with a given tolerance
tol = 1e-8;
res = abs(x1-x2)<tol;


