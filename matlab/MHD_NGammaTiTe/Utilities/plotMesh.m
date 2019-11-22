function varargout = plotMesh(X,T,eltype,option,nodesNum)

% plotMesh(X,T,faceNodes,option,nodesNum) plots the mesh defined by X and T
%
% Input:
%   X: nodal coordinates
%   T: connectivities (elements)
%   faceNodes: nOfFaces x nOfNodesPerFace matrix. Indicates the order of
%              the face nodes in a reference element. The numbering of
%              the faces has to be given in a clockwise sense, for instance:
%
%                       QUA elements                 TRI elements
%
%                             3
%                        *-------*                     *
%                         |       |                      / \
%                      4 |       | 2               3 /   \ 2
%                         |       |                    /     \
%                         *-------*                *-------*
%                             1                           1
%
%              For a given face, the column index of the matrix indicates
%              the global position of the node. This global numbering has
%              to ascend in a clockwise sense too.
%
%   option (optional): type 'plotNodes' to see the nodes' position on the
%                      ploted mesh, or type 'plotNodesNum' to see their global
%                      number.
%   nodesNum (necesary if option = 'plotNodesNum'): type 'all' to plot the
%                                                   global postion of all
%                                                   nodes, or enter a list
%                                                   array with the selected
%                                                   nodes.
%
% Output:
%   patchHandle (optional): handle to the created patch object

if eltype==0
    faceNodes = faceNodes_aux_quads(size(T,2));
    optn = 1;
else
    faceNodes = faceNodes_aux(size(T,2));
    optn = 1;
end
%Ordering the face nodes in a row vector without connectivity between them
[nOfFaces,nOfNodesPerFace] = size(faceNodes);
oFaceNodes = zeros(1,nOfFaces*(nOfNodesPerFace-1));
np = nOfNodesPerFace - 1;
aux = 1 - np;
aux2 = 0;
for iface = 1:nOfFaces
    aux = aux + np;
    aux2 = aux2 + np;
    oFaceNodes(aux:aux2) = faceNodes(iface,1:np);
end

%Conectivity for the faces
patchFaces = T(:,oFaceNodes);

%Plot mesh
patchHandle = patch('Faces',patchFaces,'Vertices',X,'FaceColor','none','EdgeAlpha',1,'linewidth',0.4);
axis equal



%Optional plots
if nargin > (2+optn) && ischar(option)
    hold on
    if strcmpi(option,'plotNodes')
        plot(X(:,1),X(:,2),'o','markerSize',3,'markerFaceColor','b')
    elseif (nargin == (4+optn)) && strcmpi(option,'plotNodesNum')
        if strcmpi(nodesNum,'all')
            list = 1:size(X,1);
            fontSize = 10;
        elseif ~isnumeric(nodesNum)
            error('wrong list of nodes')
        else
            list = nodesNum;
            fontSize = 15;
            plot(X(list,1),X(list,2),'o','markerSize',3,'markerFaceColor','b')
        end
        for inode = list
            text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize,...
                'Color',[1 0 0])
        end
    elseif (nargin == (4+optn)) && strcmpi(option,'plotNodesNumAndElements')
        if strcmpi(nodesNum,'all')
            list = 1:size(X,1);
            fontSize = 16;
        elseif ~isnumeric(nodesNum)
            error('wrong list of nodes')
        else
            list = nodesNum;
            fontSize = 15;
            plot(X(list,1),X(list,2),'o','markerSize',3,'markerFaceColor','b')
        end
        for inode = list
            text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize,...
                'Color',[1 0 0])
        end
        for iElem = 1:size(T,1)
            xbar = 1/3*(X(T(iElem,1),1)+X(T(iElem,2),1)+X(T(iElem,3),1));
            ybar = 1/3*(X(T(iElem,1),2)+X(T(iElem,2),2)+X(T(iElem,3),2));
            text(xbar,ybar,int2str(iElem),'FontSize',fontSize+2,...
                'Color',[0 0 1])
        end

    elseif (nargin == (3+optn)) && strcmpi(option,'plotElements')
        fontSize = 15;
        for iElem = 1:size(T,1)
            xbar = 1/3*(X(T(iElem,1),1)+X(T(iElem,2),1)+X(T(iElem,3),1));
            ybar = 1/3*(X(T(iElem,1),2)+X(T(iElem,2),2)+X(T(iElem,3),2));
            text(xbar,ybar,int2str(iElem),'FontSize',fontSize+2,...
                'Color',[0 0 1])
        end  
    else
        error('wrong optional argument. Check help to fix the error')
    end
    hold off
end

%Output variable
if ~nargout
    varargout = {};
else
    varargout = {patchHandle};
end

function res = faceNodes_aux(nOfElementNodes)
switch nOfElementNodes
    case 3 %P1
        res = [1 2; 2 3; 3 1];
    case 6 %P2
        res = [1 4 2; 2 5 3; 3 6 1];
    case 10 %P3
        res = [1 4 5 2; 2 6 7 3; 3 8 9 1];
    case 15 %P4
        res = [1 4 5 6 2; 2 7 8 9 3; 3 10 11 12 1];
    case 21 %P5
        res = [1 4:7 2; 2 8:11 3; 3 12:15 1];
    case 28 %P6
        res = [1 4:8 2; 2 9:13 3; 3 14:18 1];
    case 36 %P7
        res = [1 4:9 2; 2 10:15 3; 3 16:21 1];
    case 45 %P8
        res = [1 4:10 2; 2 11:17 3; 3 18:24 1];
    case 55 %P9
        res = [1 4:11 2; 2 12:19 3; 3 20:27 1];
    case 66 %P10
        res = [1 4:12 2; 2 13:21 3; 3 22:30 1];
    case 78 %P11
        res = [1 4:13 2; 2 14:23 3; 3 24:33 1];
end

 function res = faceNodes_aux_quads(nOfElementNodes)
switch nOfElementNodes
    case 4 %P1
        res = [1 2; 2 3; 3 4; 4 1];
    case 9 %P2
        res = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
    case 16 %P3
        res = [1 5 6 2; 2 7 8 3; 3 9 10 4; 4 11 12 1]; 
   case 25 %P4
        res = [1 5 6  7 2; 2 8 9 10 3; 3 11 12 13 4; 4 14 15 16 1];    
   case 36 %P5
        res = [1 5 6  7  8 2; 2 9 10 11 12 3; 3 13 14 15 16 4; 4 17 18 19 20 1];            
   case 49 %P6
        res = [1 5 6  7  8 9 2; 2 10 11 12 13 14 3; 3 15 16 17 18 19  4; 4 20 21 22 23 24  1];            
   case 64 %P7
        res = [1 5 6  7  8 9 10 2; 2 11 12 13 14 15 16 3; 3 17 18 19 20 21 22 4; 4 23 24 25 26 27 28 1];            
   case 81 %P8
        res = [1 5 6  7  8 9 10 11 2; 2 12 13 14 15 16 17 18 3; 3 19 20 21 22 23 24 25 4; 4  26 27 28 29 30 31 32 1];            
        
end


