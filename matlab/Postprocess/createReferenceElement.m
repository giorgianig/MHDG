function theReferenceElement=createReferenceElement(elementType,nOfElementNodes,handle,varargin)

% theReferenceElement=createReferenceElement(elementType,nOfElementNodes)
% Input:
%  elementType: 0 for quadrilateral, 1 for triangle
%  nOfElementNodes: number of nodes in the reference element
%  nOfGaussPoints (optional): n� of gauss points of the 1D quadrature. The
%                             default value is nDeg + 2 where nDeg is the
%                             degree of interpolation
% Output:
%  theReferenceElement: struct containing
%     .IPcoordinates: coordinates of the integration points for 2D elemens
%     .IPweights: weights of the integration points for 2D elements
%     .N: shape functions at the IP
%     .Nxi,.Neta: derivatives of the shape functions at the IP
%     .IPcoordinates1d: coordinates of the integration points for 1D boundary elemens
%     .IPweights1d: weights of the integration points for 1D boundary elements
%     .N1d: 1D shape functions at the IP
%     .N1dxi: derivatives of the 1D shape functions at the IP
%     .faceNodes: matrix [nOfFaces nOfNodesPerFace] with the edge nodes numbering
%     .innerNodes: vector [1 nOfInnerNodes] with the inner nodes numbering
%     .faceNodes1d: vector [1 nOfNodesPerElement] with the 1D nodes numbering
%     .NodesCoord: spatial coordinates of the element nodes
%     .NodesCoord1d: spatial coordinates of the 1D element nodes

if elementType == 1 %Triangles
    
    switch nOfElementNodes
        case 3 %P1
            nDeg = 1;
            faceNodes = [1 2; 2 3; 3 1];
            innerNodes = [];
            faceNodes1d = 1:2;
            coord2d = [-1 -1; 1 -1; -1 1];
            coord1d = [-1; 1];
        case 6 %P2
            nDeg = 2;
            faceNodes = [1 4 2; 2 5 3; 3 6 1];
            innerNodes = [];
            faceNodes1d = 1:3;
            coord2d = [-1 -1; 1 -1; -1 1; 0 -1; 0 0; -1 0];
            coord1d = [-1; 0; 1];
        case 10 %P3
            nDeg = 3;
            faceNodes = [1 4 5 2; 2 6 7 3; 3 8 9 1];
            innerNodes = 10;
            faceNodes1d = 1:4;
        case 15 %P4
            nDeg = 4;
            faceNodes = [1 4 5 6 2; 2 7 8 9 3; 3 10 11 12 1];
            innerNodes = 13:15;
            faceNodes1d = 1:5;
        case 21 %P5
            nDeg = 5;
            faceNodes = [1 4:7 2; 2 8:11 3; 3 12:15 1];
            innerNodes = 16:21;
            faceNodes1d = 1:6;
        case 28 %P6
            nDeg = 6;
            faceNodes = [1 4:8 2; 2 9:13 3; 3 14:18 1];
            innerNodes = 19:28;
            faceNodes1d = 1:7;
        case 36 %P7
            nDeg = 7;
            faceNodes = [1 4:9 2; 2 10:15 3; 3 16:21 1];
            innerNodes = 22:36;
            faceNodes1d = 1:8;
        case 45 %P8
            nDeg = 8;
            faceNodes = [1 4:10 2; 2 11:17 3; 3 18:24 1];
            innerNodes = 25:45;
            faceNodes1d = 1:9;
        case 55 %P9
            nDeg = 9;
            faceNodes = [1 4:11 2; 2 12:19 3; 3 20:27 1];
            innerNodes = 28:55;
            faceNodes1d = 1:10;
        case 66 %P10
            nDeg = 10;
            faceNodes = [1 4:12 2; 2 13:21 3; 3 22:30 1];
            innerNodes = 31:66;
            faceNodes1d = 1:11;
        case 78 %P11
            nDeg = 11;
            faceNodes = [1 4:13 2; 2 14:23 3; 3 24:33 1];
            innerNodes = 34:78;
            faceNodes1d = 1:12;
        otherwise
            error('Not implemented yet')
    end
    
    if nDeg >= 3 %EZ4U rules imply fekete nodes
        %coord2d = feketeNodesTri2D(nDeg,faceNodes,innerNodes);
        feketeFile = load('positionFeketeNodesTri2D_EZ4U.mat'); %EZ4U reference element
        coord2d = feketeFile.feketeNodesPosition.(['P' num2str(nDeg)]);
        coord1d = feketeNodes1D(nDeg,faceNodes1d);
%                         [coord2d coord1d] = equispacedNodesReferenceElement(nDeg);  % for equispaced nodes.. (20/12/2010 G.G.)
    end
    
elseif elementType == 0 %Quadrilaterals
    
    switch nOfElementNodes
        case 4 %Q1
            nDeg = 1;
            faceNodes = [1 2; 2 3; 3 4; 4 1];
            innerNodes = [];
            faceNodes1d = 1:2;
            coord2d = [-1 -1; 1 -1; 1 1; -1 1];
            coord1d = [-1; 1];
            perm = [1 2 4 3];
        case 9 %Q2
            nDeg = 2;
            faceNodes = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
            innerNodes = 9;
            faceNodes1d = 1:3;
            coord2d = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];
            coord1d = [-1; 0; 1];
            perm = [ 1     5     2     8     9     6     4     7     3];
        case 16 %Q3
            nDeg = 3;
            faceNodes =[1 5 6 2; 2 7 8 3; 3 9 10 4; 4 11 12 1]; 
            innerNodes = 13:16;
            faceNodes1d = 1:4;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [ 1     5     6     2    12    13    14     7    11    15    16     8     4    10     9     3];
            coord2d(perm,:) = coord2d;
 
         case 25 %Q4
            nDeg = 4;
            faceNodes =[1 5 6  7 2; 2 8 9 10 3; 3 11 12 13 4; 4 14 15 16 1];    
            innerNodes = 17:25;
            faceNodes1d = 1:5;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [ 1     5     6     7     2    16    17    18    19     8    15    20    21  22     9    14    23    24    25    10     4    13    12    11     3 ];
            coord2d(perm,:) = coord2d;
         case 36 %Q5
            nDeg = 5;
            faceNodes = [1 5 6  7  8 2; 2 9 10 11 12 3; 3 13 14 15 16 4; 4 17 18 19 20 1];      
            innerNodes = 21:36;
            faceNodes1d = 1:6;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [     1     5     6     7     8     2    20    21    22    23    24     9    19  25    26    27    28    10    18    29    30    31    32    11    17    33  34    35    36    12     4    16    15    14    13     3 ];
            coord2d(perm,:) = coord2d;            
         case 49 %Q6
            nDeg = 6;
            faceNodes = [1 5 6  7  8 9 2; 2 10 11 12 13 14 3; 3 15 16 17 18 19  4; 4 20 21 22 23 24  1];     
            innerNodes = 25:49;
            faceNodes1d = 1:7;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [      1     5     6     7     8     9     2    24    25    26    27    28    29   10    23    30    31    32    33    34    11    22    35    36    37    38 39    12    21    40    41    42    43    44    13    20    45    46    47 48    49    14     4    19    18    17    16    15     3];
            coord2d(perm,:) = coord2d;            

         case 64 %Q7
            nDeg = 7;
            faceNodes = [1 5 6  7  8 9 10 2; 2 11 12 13 14 15 16 3; 3 17 18 19 20 21 22 4; 4 23 24 25 26 27 28 1];           
            innerNodes = 29:64;
            faceNodes1d = 1:8;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [  1     5     6     7     8     9    10     2    28    29    30    31    32 33    34    11    27    35    36    37    38    39    40    12    26    41 42    43    44    45    46    13    25    47    48    49    50    51    52  14    24    53    54    55    56    57    58    15    23    59    60    61   62    63    64    16     4    22    21    20    19    18    17     3];
            coord2d(perm,:) = coord2d;            

          case 81 %Q8
            nDeg = 8;
            faceNodes =  [1 5 6  7  8 9 10 11 2; 2 12 13 14 15 16 17 18 3; 3 19 20 21 22 23 24 25 4; 4  26 27 28 29 30 31 32 1];             
            innerNodes = 33:81;
            faceNodes1d = 1:9;
            coord1d = feketeNodes1D(nDeg,faceNodes1d);
            coord2d = zeros((nDeg+1)^2,2);
            for i = 1:nDeg+1
                ind = (i-1)*(nDeg+1)+(1:(nDeg+1));
                coord2d(ind,1) = coord1d;
                coord2d(ind,2) = coord1d(i);
            end
            perm = [   1     5     6     7     8     9    10    11     2    32    33    34    35 36    37    38    39    12    31    40    41    42    43    44    45    46  13    30    47    48    49    50    51    52    53    14    29    54    55  56    57    58    59    60    15    28    61    62    63    64    65    66  67    16    27    68    69    70    71    72    73    74    17    26    75  76    77    78    79    80    81    18     4    25    24    23    22    21  20    19     3];
            coord2d(perm,:) = coord2d;            
        otherwise
            error('Not implemented yet')
    end
    
else
    error('Element not allowed')
end

%Compute shape functions and quadrature
if isempty(varargin)
        nOfGaussPoints = nDeg + 2;
%     nOfGaussPoints = 2*nDeg + 1;
else
    nOfGaussPoints = varargin{:};
end
if elementType == 1
    [shapeFun2d,gw2d,gp2d] = ...
        computeShapeFunctionsReferenceElement(nDeg,coord2d,nOfGaussPoints,elementType);
    [shapeFun1d,gw1d,gp1d] = ...
        computeShapeFunctionsReferenceElement(nDeg,coord1d,nOfGaussPoints);
elseif elementType == 0
    [shapeFun1d,gw1d,gp1d] = ...
        computeShapeFunctionsReferenceElement(nDeg,coord1d,nOfGaussPoints);
    [shapeFun2d,gw2d,gp2d] = createShapeFunctions2dTensor(shapeFun1d,gw1d,gp1d,perm);
%         [shapeFun2d,gw2d,gp2d] = ProvaComputeShapeFunctionsReferenceElementQuads(nDeg,coord1d,nOfGaussPoints);

end
N = shapeFun2d(:,:,1)';
Nxi = shapeFun2d(:,:,2)';
Neta = shapeFun2d(:,:,3)';
N1 = shapeFun1d(:,:,1)';
Nxi1 = shapeFun1d(:,:,2)';

%Creating reference element structure
theReferenceElement = struct('IPcoordinates',gp2d,...
    'IPweights',gw2d,'N',N,'Nxi',Nxi,'Neta',Neta,...
    'IPcoordinates1d',gp1d,'IPweights1d',gw1d,...
    'N1d',N1,'N1dxi',Nxi1,'faceNodes',faceNodes,...
    'innerNodes',innerNodes,'faceNodes1d',faceNodes1d,...
    'NodesCoord',coord2d,'NodesCoord1d',coord1d,'degree',nDeg,...
    'elemType',elementType);

