% function GenerateHDF5_meshes_forParallel()

for i = 1:ndiv
    if strcmpi(meshName(end-2:end),'.h5')
        fileName = [meshName(1:end-3) '_' num2str(i) '_' num2str(ndiv) '.h5'];
    else
        fileName = [meshName(1:end-4) '_' num2str(i) '_' num2str(ndiv) '.h5'];
    end 
    X = XX{i};
    T = TT{i};
    loc2glob_el = Loc2Glob_el{i};
    loc2glob_fa = Loc2Glob_fa{i};
    loc2glob_no = Loc2Glob_no{i};
    ghostFaces  = GhostFaces{i};
    ghostElems  = GhostElems{i};
    flipFace    = flipFaces_GHOST{i};
    ghostPro    = GhostPro{i};
    ghostLoc    = GhostLoc{i};
    ghelsPro    = GhelsPro{i};
    ghelsLoc    = GhelsLoc{i};
    nghels      = sum(ghostElems);
    nf = size(loc2glob_fa,1);
    for ibound = 1:nboundaries
        name = boundaryNames{ibound};
        eval(['Tb_' name(4:end) ' = TTb_' name(4:end) '{i};'])
    end
    Tb_CON = TTb_CON{i};
    
    loc2glob = Loc2Glob_no{i};
    
    % Starting conversion
    delete(fileName)
    Nnodes = size(X,1);
    Nelems = size(T,1);
    Ndim   = size(X,2);
    Nnodesperelem = size(T,2);
    h5create(fileName,'/T',size(T));
    h5create(fileName,'/X',size(X));
    h5create(fileName,'/Nnodes',1);
    h5create(fileName,'/Ndim',1);
    h5create(fileName,'/Nelems',1);
    h5create(fileName,'/elemType',1);
    h5create(fileName,'/Nnodesperelem',1);
    % Boundaries
    s1 = 0;
    boundInt = zeros(nboundaries,1);
    for ibound = 1:nboundaries
        name = boundaryNames{ibound};
        eval(['aux_bound =  TTb_' name(4:end) '{i};'])
        s1 = s1+size(aux_bound,1);
        if any(strcmpi(name,{'Tb_Dirichlet','Tb_Diriclet'}))
            boundInt(ibound) = 1;
        elseif strcmpi(name,'Tb_LEFT')
            boundInt(ibound) = 2;
        elseif strcmpi(name,'Tb_RIGHT')
            boundInt(ibound) = 3;
        elseif strcmpi(name,'Tb_UP')
            boundInt(ibound) = 4;
        elseif strcmpi(name,'Tb_DOWN')
            boundInt(ibound) = 5;
        elseif strcmpi(name,'Tb_WALL')
            boundInt(ibound) = 6;
        elseif strcmpi(name,'Tb_LIM')
            boundInt(ibound) = 7;
        elseif strcmpi(name,'Tb_IN')
            boundInt(ibound) = 8;
        elseif strcmpi(name,'Tb_OUT')
            boundInt(ibound) = 9;
        elseif strcmpi(name,'Tb_ULIM')
            boundInt(ibound) = 10;
        end
    end
    aux_bound = Tb_CON;
    s1 = s1+size(aux_bound,1);
    s2 = size(aux_bound,2);
    Tb = zeros(s1,s2);
    boundaryFlag = zeros(s1,1);
    s1 = 1;
    for ibound = 1:nboundaries
        name = boundaryNames{ibound};
        eval(['aux_bound = TTb_' name(4:end) '{i};'])
        Tb(s1:size(aux_bound,1)+s1-1,:) = aux_bound;
        boundaryFlag(s1:size(aux_bound,1)+s1-1) = boundInt(ibound);
        s1=s1+size(aux_bound,1);
    end
    aux_bound = Tb_CON;
    Tb(s1:size(aux_bound,1)+s1-1,:) = aux_bound;
    h5create(fileName,'/Nnodesperface',1);
    h5create(fileName,'/Nextfaces',1);
    h5create(fileName,'/Nfaces',1);
    h5create(fileName,'/Tb',size(Tb));
    h5create(fileName,'/boundaryFlag',size(boundaryFlag));
    h5create(fileName,'/loc2glob_fa',nf);
    h5create(fileName,'/loc2glob_el',Nelems);
    h5create(fileName,'/loc2glob_no',Nnodes);
    h5create(fileName,'/ghostFaces',nf);
    h5create(fileName,'/ghostElems',size(T,1));
    h5create(fileName,'/ghostFlp',size(flipFace));
    h5create(fileName,'/ghostPro',size(flipFace))
    h5create(fileName,'/ghostLoc',size(flipFace))
    if nghels==0
        h5create(fileName,'/ghelsPro',1)
        h5create(fileName,'/ghelsLoc',1)
    else
        h5create(fileName,'/ghelsPro',nghels)
        h5create(fileName,'/ghelsLoc',nghels)
    end
    
    % Write in the database
    h5disp(fileName)
    h5write(fileName,'/T',T);
    h5write(fileName,'/X',X)
    h5write(fileName,'/Nnodes',Nnodes)
    h5write(fileName,'/Ndim',Ndim)
    h5write(fileName,'/Nelems',Nelems)
    h5write(fileName,'/Nfaces',nf)
    h5write(fileName,'/Nnodesperelem',Nnodesperelem)
    h5write(fileName,'/elemType',elemTypeFort)
    h5write(fileName,'/Nextfaces',size(Tb,1))
    h5write(fileName,'/Nnodesperface',s2)
    h5write(fileName,'/Tb',Tb)
    h5write(fileName,'/boundaryFlag',boundaryFlag)
    h5write(fileName,'/loc2glob_fa',loc2glob_fa)
    h5write(fileName,'/loc2glob_el',loc2glob_el)
    h5write(fileName,'/loc2glob_no',loc2glob_no)
    h5write(fileName,'/ghostFaces',ghostFaces)
    h5write(fileName,'/ghostElems',ghostElems)
    h5write(fileName,'/ghostFlp',flipFace)
    h5write(fileName,'/ghostPro',ghostPro)
    h5write(fileName,'/ghostLoc',ghostLoc)
    if nghels==0
        h5write(fileName,'/ghelsPro',-1)
        h5write(fileName,'/ghelsLoc',-1)
    else
        h5write(fileName,'/ghelsPro',ghelsPro)
        h5write(fileName,'/ghelsLoc',ghelsLoc)
    end
    % Check
    h5disp(fileName)
    
end