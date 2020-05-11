function varargout = GenerateMatFileFromEZ4U(arg1,arg2,arg3,arg4,arg5)

% This is the same function that GenerateMatFileFromEZ4U but the boundary
% connectivity matrices. These incloude the element and the face.

%Checking output variable's call
if nargout > 1
    error(['You are calling more than one output variable. There is only '...
        'a cell array as an output variable'])
end


%Checking input syntax and detecting EZ4U mesh files
switch nargin
    case 1
        if isCellStringOrChar(arg1)
            meshFiles = getEZ4UMeshFiles(arg1);
        elseif isNumVector(arg1)
            meshFiles = getEZ4UMeshFiles();
            face_nodes1d = arg1;
        else error('Wrong input arguments. Check help to correctly use the syntax')
        end
    case 2
        if isCellStringOrChar(arg1,arg2)
            meshFiles = getEZ4UMeshFiles(arg1,arg2);
        elseif isNumMatrix(arg1) && isVector(arg2)
            meshFiles = getEZ4UMeshFiles();
            face_nodes = arg1; inner_nodes = arg2;
        elseif isCellStringOrChar(arg1) && isNumVector(arg2)
            meshFiles = getEZ4UMeshFiles(arg1);
            face_nodes1d = arg2;
        else error('Wrong input arguments. Check help to correctly use the syntax')
        end
    case 3
        if isCellStringOrChar(arg1,arg2) && isNumVector(arg3)
            meshFiles = getEZ4UMeshFiles(arg1,arg2);
            face_nodes1d = arg3;
        elseif isNumMatrix(arg1) && isVector(arg2) && isNumVector(arg3)
            meshFiles = getEZ4UMeshFiles();
            face_nodes = arg1; inner_nodes = arg2; face_nodes1d = arg3;
        elseif isCellStringOrChar(arg1) && isNumMatrix(arg2) && isVector(arg3)
            meshFiles = getEZ4UMeshFiles(arg1);
            face_nodes = arg2; inner_nodes = arg3;
        else error('Wrong input arguments. Check help to correctly use the syntax')
        end
    case 4
        if isCellStringOrChar(arg1,arg2) && isNumMatrix(arg3) && isVector(arg4)
            meshFiles = getEZ4UMeshFiles(arg1,arg2);
            face_nodes = arg3; inner_nodes = arg4;
        elseif isCellStringOrChar(arg1) && isNumMatrix(arg2) && isVector(arg3)...
                && isNumVector(arg4)
            meshFiles = getEZ4UMeshFiles(arg1);
            face_nodes = arg2; inner_nodes = arg3; face_nodes1d = arg4;
        else error('Wrong input arguments. Check help to correctly use the syntax')
        end
    case 5
        if isCellStringOrChar(arg1,arg2) && isNumMatrix(arg3) && isVector(arg4)...
                && isNumVector(arg5)
            meshFiles = getEZ4UMeshFiles(arg1,arg2);
            face_nodes = arg3; inner_nodes = arg4; face_nodes1d = arg5;
        else error('Wrong input arguments. Check help to correctly use the syntax')
        end
    otherwise
        meshFiles = getEZ4UMeshFiles();
end
nOfMeshFiles = numel(meshFiles);


%Reading mesh files and creating .mat files
matFiles = cell(1,nOfMeshFiles);
for ifile = 1:nOfMeshFiles
    ez4uFile = meshFiles{ifile};

    %Open EZ4U mesh file
    fid = fopen(ez4uFile,'r');

    %Reading header
    header = fscanf(fid,'%d',3);
    nOfNodes = header(1);
    nOfElements = header(2);
    nOfAttributes = header(3);

    %Reading nodes
    dim = 2;
    X = zeros(nOfNodes,dim);
    for iNode = 1:nOfNodes
        lineVar = fscanf(fid,'%f',4);
        X(iNode,:) = lineVar(2:dim+1);
    end

    %Element information
    nOfInfo = 7;
    Info = fscanf(fid,'%d',nOfInfo);
    nOfFaces = Info(4);
    nOfEdgeNodes = Info(5);
    nOfInnerNodes = Info(6);
    nOfElementNodes = nOfInnerNodes + nOfEdgeNodes*nOfFaces + nOfFaces;
    elemType = 4 - nOfFaces;

    %Reading connectivities for all the mesh
    T = zeros(nOfElements,nOfElementNodes);
    repeatedFaceInfo = zeros(nOfElements,2);
    firstRow = fscanf(fid,'%d',nOfElementNodes);
    T(1,:) = firstRow';
    lineSize = nOfInfo + nOfElementNodes;
    nOfElemNodesPos = nOfInfo+1:nOfElementNodes+nOfInfo;
    for iElem = 2:nOfElements
        lineVar = fscanf(fid,'%d',lineSize);
        T(iElem,:) = lineVar(nOfElemNodesPos)';
    end

    %Reading connectivities for the boundary and submeshes (attributes)
    indexSubMesh = 0;
    indexBoundary = 0;
    newIndex = nOfElements;
    if nOfAttributes
        attribInfo = cell(2,nOfAttributes,3);
        for iAttrib = 1:nOfAttributes
            lineIntegers = fscanf(fid,'%d',4);
            nOfMarkedElements = lineIntegers(4);
            attribType = lineIntegers(2);
            attributeName = fscanf(fid,'%s',1);

            if attribType == 2 %Submeshes (2D faces)
                indexSubMesh = indexSubMesh + 1;
                attribInfo{1,indexSubMesh,1} = zeros(nOfMarkedElements,nOfElementNodes);
                for item = 1:nOfMarkedElements
                    lineVar = fscanf(fid,'%d',3);
                    element = lineVar(1);
                    attribInfo{1,indexSubMesh,1}(item,:) = T(element,:);
                end

                %Attribute's name
                attribInfo{1,indexSubMesh,2} = ['T_' attributeName];

            elseif attribType == 1 %Boundaries (2D edges)
                indexBoundary = indexBoundary + 1;
                attribInfo{2,indexBoundary,1} = zeros(nOfMarkedElements,nOfEdgeNodes+2);
                attribInfo{2,indexBoundary,3} = zeros(nOfMarkedElements,2);
                for item = 1:nOfMarkedElements
                    lineVar = fscanf(fid,'%d',3);
                    element = lineVar(1);
                    face = lineVar(2);
                    elementNodes = [T(element,1:nOfFaces) zeros(1,nOfEdgeNodes)...
                        T(element,nOfFaces+1:end)];
                    elementFaces = [elementNodes(1:nOfFaces) elementNodes(1)];
                    vertexNodes = elementFaces(face:face+1);
                    ini = nOfFaces + face*nOfEdgeNodes + 1;
                    fin = ini + nOfEdgeNodes - 1;
                    innerNodes = elementNodes(ini:fin);
                    attribInfo{2,indexBoundary,1}(item,:) = [vertexNodes(1) innerNodes vertexNodes(2)];
                    if repeatedFaceInfo(element,1)
                        oldIndex = repeatedFaceInfo(element,2);
                        attribInfo{2,indexBoundary,3}(item,:) = [oldIndex face];
                    else
                        repeatedFaceInfo(element,:) = [1 newIndex];
                        attribInfo{2,indexBoundary,3}(item,:) = [newIndex face];
                        newIndex = newIndex - 1;
                    end
                end

                %Attribute's name
                attribInfo{2,indexBoundary,2} = ['Tb_' attributeName];

            else
                error('Attribute type not implemented')
            end
        end
    end

    %Renumerate elements in order to put interior elements first
    extElem = find(repeatedFaceInfo(:,1) == 1);
    permutations = repeatedFaceInfo(extElem,2);
    [aux1,aux2] = sort(permutations);
    extElem = extElem(aux2);
    intElem = setdiff((1:nOfElements)',extElem);
    T = T([intElem;extElem],:);

    %Applying the given local numbering to the boundary (if it's necessary)
    if exist('face_nodes1d','var')
        nOfNodes1d = size(attribInfo{2,1,1},2);
        givenNOfNodes1d = length(face_nodes1d);
        if (nOfNodes1d ~= givenNOfNodes1d) && nOfNodes1d
            error(['Wrong argument face_nodes1d for the file "' ez4uFile...
                '". You specified a wrong number of nodes'])
        end
        for item = 1:indexBoundary
            attribInfo{2,item,1}(:,face_nodes1d) = attribInfo{2,item,1};
        end

        %info to store in elemInfo
        elemFaceNodes1d = face_nodes1d;
    else
        elemFaceNodes1d = 1:nOfEdgeNodes+2;
    end

    %Applying the given local numbering to the mesh and submeshes (if it's necessary)
    if exist('face_nodes','var') && exist('inner_nodes','var')
        [givenNOfFaces,givenNOfFaceNodes] = size(face_nodes);
        givenNOfInnerNodes = length(inner_nodes);
        if givenNOfFaceNodes ~= nOfEdgeNodes+2
            error(['Wrong argument face_nodes for the file "' ez4uFile...
                '". You specified a wrong number of face nodes'])
        elseif givenNOfFaces ~= nOfFaces
            error(['Wrong argument face_nodes for the file "' ez4uFile...
                '". You specified a wrong number of faces'])
        elseif givenNOfInnerNodes ~= nOfInnerNodes
            error(['Wrong argument inner_nodes for the file "' ez4uFile...
                '". You specified a wrong number of inner nodes'])
        end
        givenOrder = getEZ4UOrder(face_nodes,inner_nodes,nOfFaces,nOfEdgeNodes);

        %Complete mesh
        T(:,givenOrder) = T;

        %Submeshes
        for item = 1:indexSubMesh
            attribInfo{1,item,1}(:,givenOrder) = attribInfo{1,item,1};
        end

        %info to store in elemInfo
        elemFaceNodes = face_nodes;
    else
        elemFaceNodes = zeros(nOfFaces,nOfEdgeNodes+2);
        elemVertices = 1:nOfFaces;
        elemEdgeNodes = nOfFaces+1:nOfEdgeNodes*nOfFaces + nOfFaces;
        elemFaceNodes(:,1) = elemVertices';
        elemFaceNodes(:,end) = [elemVertices(2:end) 1]';
        aux = 0;
        for icolumn = 2:nOfEdgeNodes+1
            ini = 1 + aux;
            fin = ini + nOfEdgeNodes*(nOfFaces - 1);
            elemFaceNodes(:,icolumn) = elemEdgeNodes(ini:nOfEdgeNodes:fin)';
            aux = aux + 1;
        end
    end

    %Creating submeshes and boundary conectivities in the current workspace
    index = [indexSubMesh indexBoundary];
    elementFaceInfo = struct();
    for iattrib = 1:length(index)
        for item = 1:index(iattrib)
            name = attribInfo{iattrib,item,2};
            var2save = genvarname(name);
            evalc([var2save '=attribInfo{iattrib,item,1}']);
            if iattrib == 2
                elementFaceInfo.(name(4:end)) = attribInfo{iattrib,item,3};
            end
        end
    end

    %Store the element information into a structure
    elemInfo = struct('type',elemType,'nOfNodes',nOfElementNodes,...
        'faceNodes',elemFaceNodes,'faceNodes1d',elemFaceNodes1d);

    %Close EZ4U mesh file
    fclose(fid);

    %Creating .mat file
    matFileName = [ez4uFile(1:end-4) '.mat'];
    matFiles{ifile} = matFileName;
    save(matFileName,'X','T*','elemInfo','elementFaceInfo')
    clear X T* elemInfo elementFaceInfo
end


%Assign output variable
if nargout == 1, varargout = {matFiles};
else varargout = [];
end


%_______________________________________________________________________
% HELPER FUNCTIONS TO CHECK THE KIND OF ARRAY THE DATA IS AND MAKE THIS
% CHECKING MORE READABLE

function isIt = isCellStringOrChar(varargin)

isIt = true;
for i = 1:length(varargin)
    data = varargin{i};
    if iscell(data)
        for j = 1:numel(data)
            if ~ischar(data{j}), isIt = false; return, end
        end
    elseif ischar(data)
    else isIt = false; return
    end
end

%---------

function isIt = isNumMatrix(data)

if isnumeric(data) && all(size(data) > 1), isIt = true;
else isIt = false;
end

%---------

function isIt = isVector(data)

if isnumeric(data)
    condition1 = size(data) == 1;
    condition2 = ~size(data);
    if any(condition1 | condition2), isIt = true;
    else isIt = false;
    end
else isIt = false;
end

%---------

function isIt = isNumVector(data)

if isnumeric(data)
    condition1 = size(data) == 1;
    condition2 = size(data) > 1;
    if any(condition1) && any(condition2), isIt = true;
    else isIt = false;
    end
else isIt = false;
end


%_______________________________________________________________________
% HELPER FUNCTIONS TO FIND OUT AND/OR CHECK THE MESH FILES

function meshFiles = getEZ4UMeshFiles(arg1,arg2)

%Assign folder or files
if nargin == 1
    if ~ischar(arg1)
        checkEZ4UMeshFiles([],arg1), meshFiles = arg1; return
    else
        condition = length(arg1) > 4;
        if condition && strcmpi(arg1(end:-1:end-3),'mcd.')
            arg1 = {arg1};
            checkEZ4UMeshFiles([],arg1), meshFiles = arg1; return
        else folder = arg1;
        end
    end
elseif nargin == 2
    if ~iscell(arg2), arg2 = {arg2}; end
    checkEZ4UMeshFiles(arg1,arg2), meshFiles = arg2; return
else folder = [];
end

%This part executes the case of being only a folder
emptyFolder = isempty(folder);
if ~emptyFolder && ~exist(folder,'dir')
    error(['The specified directory "' folder '" does not exist'])
elseif ~emptyFolder, addpath(folder), folder = [folder '/'];
end
dirMeshFiles = dir([pwd '/' folder '*dcm']);
nOfDirMeshFiles = length(dirMeshFiles);
if ~nOfDirMeshFiles && ~emptyFolder
    error('There is no .dcm file into the specified directory')
elseif ~nOfDirMeshFiles && emptyFolder
    error(['There is no .dcm file into the current directory.'...
        ' You have to specify the directory where they are located'])
end
meshFiles = cell(1,nOfDirMeshFiles);
for i = 1:nOfDirMeshFiles
    meshFiles{i} = dirMeshFiles(i).name;
end

%---------

function checkEZ4UMeshFiles(folder,meshFiles)

emptyFolder = isempty(folder);
if ~emptyFolder && ~exist(folder,'dir')
    error(['The specified directory "' folder '" does not exist'])
elseif ~emptyFolder, addpath(folder)
end
nOfmeshFiles = numel(meshFiles);
for ifile = 1:nOfmeshFiles
    file = meshFiles{ifile};
    dirFile = dir([pwd '/' folder '/' file]);
    condition = length(dirFile);
    if ~condition && ~emptyFolder
        error(['The file "' file '" is not located into the'...
            ' specified directory'])
    elseif ~condition && emptyFolder
        error(['The file "' file '" is not located into the'...
            ' current directory. You have to specify the directory'...
            ' where it is located'])
    end
end


%_______________________________________________________________________
% AUXILIAR FUNCTIONS

function Order = getEZ4UOrder(fnodes,inodes,nfaces,nenodes)

vnodes = zeros(1,nfaces);
enodes = zeros(1,nfaces*nenodes);
aux = 1;
for i = 1:nfaces
    vnodes(i) = fnodes(i,1);
    enodes(aux:aux+nenodes-1) = fnodes(i,2:end-1);
    aux = aux + nenodes;
end
if size(inodes,2) == 1, inodes = inodes'; end
Order = [vnodes enodes inodes];



