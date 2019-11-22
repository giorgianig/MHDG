function varargout = plotMesh(X,T,option,nodesNum)

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
%                         *-------*                       *
%                         |       |                      / \
%                       4 |       | 2                 3 /   \ 2
%                         |       |                    /     \
%                         *-------*                   *-------*
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

faceNodes = faceNodes_aux(size(T,2));
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
% patchHandle = patch('Faces',patchFaces,'Vertices',X,'FaceColor','none','EdgeAlpha',1,'edgecolor','r');
patchHandle = patch('Faces',patchFaces,'Vertices',X,'FaceColor','none','EdgeAlpha',0.2,'linewidth',0.2);

axis equal

%Optional plots
if nargin > 2
    hold on
    if strcmpi(option,'plotNodes')
        plot(X(:,1),X(:,2),'o','markerSize',3,'markerFaceColor','n')
    elseif (nargin == 4) && strcmpi(option,'plotNodesNum')
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
    elseif (nargin == 4) && strcmpi(option,'plotNodesNumAndElements')
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

    elseif (nargin == 3) && strcmpi(option,'plotElements')
        fontSize =16;
        for iElem = 1:size(T,1)
            xbar = 1/3*(X(T(iElem,1),1)+X(T(iElem,2),1)+X(T(iElem,3),1));
            ybar = 1/3*(X(T(iElem,1),2)+X(T(iElem,2),2)+X(T(iElem,3),2));
            text(xbar,ybar,int2str(iElem),'FontSize',fontSize,...
                'Color',[0 0 1])
        end
    elseif (nargin == 3) && (isnumeric(option) || islogical(option))
        patch('Faces',patchFaces(option,:),'Vertices',X,'FaceColor','r','EdgeAlpha',1);
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
    case 91 %P12
        res = [1 4:14 2; 2 15:25 3; 3 26:36 1];
    case 105 %P13
        res = [1 4:15 2; 2 16:27 3; 3 28:39 1];
    case 120 %P14
        res = [1 4:16 2; 2 17:29 3; 3 30:42 1];
    case 136 %P15
        res = [1 4:17 2; 2 18:31 3; 3 32:45 1];
    case 153 %P16
        res = [1 4:18 2; 2 19:33 3; 3 34:48 1];
    case 171 %P17
        res = [1 4:19 2; 2 20:35 3; 3 36:51 1];
    case 190 %P18
        res = [1 4:20 2; 2 21:37 3; 3 38:54 1];
    case 210 %P19
        res = [1 4:21 2; 2 22:39 3; 3 40:57 1];
    case 231 %P20
        res = [1 4:22 2; 2 23:41 3; 3 42:60 1];
end

                  


