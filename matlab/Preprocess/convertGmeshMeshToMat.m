close all
clear all

%filename = ['Gmsh/West_Mesh_YesHole_coarse.msh'];
filename = ['Gmsh/Circ_HoriLim_YesHole_coarse.msh'];

global P 

P = 4;

fileID = fopen(filename)

%for i=1:11
% dummy = fgets(fileID);
%end 
dummy = strtrim(string(fgets(fileID)));
while dummy ~= '$Nodes'
    dummy = strtrim(string(fgets(fileID)));
end

n_nodes = fscanf(fileID,'%d',1);
format long
X = fscanf(fileID,'%g',[4,n_nodes]);
X = X(2:3,:)';

for i=1:3
    dummy = fgets(fileID);
end

n_elements = fscanf(fileID,'%d',1);

n_Tb_IN = 0;
n_Tb_OUT = 0;
n_Tb_LIM = 0;
n_T_DOM = 0;

dummy = fgets(fileID);

for i=1:n_elements
 array = str2num(fgets(fileID));
 switch array(4)
     case 1
         n_Tb_IN = n_Tb_IN + 1;
     case 2
         n_Tb_OUT = n_Tb_OUT + 1;
     case 3
         n_Tb_LIM = n_Tb_LIM + 1;
     case 4
         n_T_DOM = n_T_DOM +1;
 end
end

fclose(fileID);

fileID = fopen(filename);

Tb_IN = zeros(n_Tb_IN,P+1);
Tb_OUT = zeros(n_Tb_OUT,P+1);
Tb_LIM = zeros(n_Tb_LIM,P+1);
T = zeros(n_T_DOM,(P+1)*(P+2)/2);
i_Tb_IN = 1;
i_Tb_OUT = 1;
i_Tb_LIM = 1;
i_T_DOM = 1;

dummy = fgets(fileID);
dummy = string(dummy);

while dummy ~= num2str(n_elements)
    dummy = strtrim(string(fgets(fileID)));
end

for i=1:n_elements
    boundary = str2num(fgets(fileID));
    switch boundary(4)
        case 1
            Tb_IN(i_Tb_IN,:) = boundary(6:end);
            i_Tb_IN = i_Tb_IN + 1;
        case 2
            Tb_OUT(i_Tb_OUT,:) = boundary(6:end);
            i_Tb_OUT = i_Tb_OUT + 1;
        case 3
            Tb_LIM(i_Tb_LIM,:) = boundary(6:end);
            i_Tb_LIM = i_Tb_LIM + 1;
        case 4
            T(i_T_DOM,:) = boundary(6:end);
            i_T_DOM = i_T_DOM + 1;
    end
end

fclose(fileID)

%Reshape Tb in linear way

Matrix_Tb_OUT = Tb_OUT;
Tb_OUT(:,2:end-1) = Matrix_Tb_OUT(:,3:end);
Tb_OUT(:,end) = Matrix_Tb_OUT(:,2);

if n_Tb_LIM~=0
    Matrix_Tb_LIM = Tb_LIM;
    Tb_LIM(:,2:end-1) = Matrix_Tb_LIM(:,3:end);
    Tb_LIM(:,end) = Matrix_Tb_LIM(:,2);
end

if n_Tb_IN ~= 0 
        Matrix_Tb_IN = Tb_IN;
        Tb_IN(:,1:2) = flip(Tb_IN(:,1:2),2);
        Tb_IN(:,3:P+1) = flip(Tb_IN(:,3:P+1),2);
        Matrix_Tb_IN = Tb_IN;
        Tb_IN(:,2:end-1) = Matrix_Tb_IN(:,3:end);
        Tb_IN(:,end) = Matrix_Tb_IN(:,2);
end 
    
if P >=5
    switch P
        case 5
            [index] = [1:16,19,17,21,20,18];
        case 6
            [index] = [1:19,22,23,20,27,28,24,26,25,21];
        case 8
            [index] = [1:25,28:31,26,39,40,43,41,32,38,45,44,33,37,42,34,36,35,27];
    end
    T = T(:,index);
end


%Generate Fekete Nodes

elemType =1;
Xmod = GenereteFeketeNodes(X,T);
X = Xmod;

if n_Tb_IN ~= 0
    if contains(filename,'ITER')
        fileName = ['ITER_YesHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
    elseif contains(filename,'West')
        fileName = ['West_YesHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
    else
        if contains(filename,'/')
            fileName = [extractAfter(extractBefore(filename,'_YesHole'),'/') '_Triangs_YesHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
        else
            fileName = [extractBefore(filename,'_YesHole') '_Triangs_YesHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
        end
    end
else
    if contains(filename,'ITER')
        fileName = ['ITER_NoHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
    elseif contains(filename,'West')
        fileName = ['West_NoHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
    else
        if contains(filename,'/')
            fileName = [extractAfter(extractBefore(filename,'_NoHole'),'/') '_Triangs_NoHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
        else
            fileName = [extractBefore(filename,'_YesHole') '_Triangs_NoHole_Nel' num2str(size(T,1)) '_P' num2str(P) '.mat'];
        end
    end
end
fileName = fullfile('../Meshes',fileName);
save(fileName,'X','Tb_IN','Tb_OUT','Tb_LIM','T')
elemType = 1;

convertMatMeshToHDF5
