function Matrix = readCSRtxtParall(Mcsr,nproc)
% Read CSR matrix and convert to Matlab sparse

dimension = 0;
I = [];
J = [];
K = [];
Loc2glob = [];
for iproc = 1:nproc
    % open file
    fid = fopen([Mcsr '_' num2str(iproc-1) '_' num2str(nproc) '.txt'],'r');
    
    % read n and nnz
    n = fscanf(fid,'%d',1);
    nonz = fscanf(fid,'%d',1);
    
    col_ind = zeros(nonz,1);
    val = col_ind;
    row_ptr = zeros(n+1,1);
    loc2glob = zeros(n,1);
    
    % skip one
    aux = fscanf(fid,'%s',1);
    % read vals
    for i = 1:nonz
        val(i) = fscanf(fid,'%f',1);
    end
    
    % skip one
    aux = fscanf(fid,'%s',1);
    % read cols
    for i = 1:nonz
        col_ind(i) = fscanf(fid,'%d',1);
    end
    
    % skip one
    aux = fscanf(fid,'%s',1);
    % read rows
    for i = 1:n+1
        row_ptr(i) = fscanf(fid,'%d',1);
    end
    
    % skip one
    aux = fscanf(fid,'%s',1);
    % read loc2glob
    for i = 1:n
        loc2glob(i) = fscanf(fid,'%d',1);
    end
    
    % momentary fix sequential
    % loc2glob = 1:n;
    
    % close file
    fclose(fid);
    
    % Conversion
    Iloc = zeros(nonz,1);
    Jloc = Iloc;
    Kloc = Iloc;
    
if any(diff(loc2glob)<0)
    disp(['Inversione nel processo: ', num2str(iproc)])
%     loc2glob
%     stop
end    
    Loc2glob = [Loc2glob;loc2glob];
        
    for r = 1:numel(row_ptr)-1
        
        for i = row_ptr(r):(row_ptr(r+1)-1)
            
            Iloc(i) = loc2glob(r);
            Jloc(i) = col_ind(i);
            Kloc(i) = val(i);
        end
    end
    
% if any(loc2glob>7536), stop, end   
% if any(col_ind>7536), stop,end
Matrix_ipro = sparse(Iloc(Iloc~=0),Jloc(Iloc~=0),Kloc(Iloc~=0));    
figure, spy(Matrix_ipro),axis([0,7536,0,7536]),axis('equal')
    
    I = [I; Iloc];
    J = [J; Jloc];
    K = [K; Kloc];
    dimension = dimension+n;   
end
Matrix = sparse(I(I~=0),J(I~=0),K(I~=0));
