clear
close all

global Mesh
kT=25;
lscale =  1.901*1e-3;
nproc = 4;

matmeshpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
solpath = '/home/giorgio/Dropbox/Fortran/MHDG_3D/test/';
solution = 'Sol2D_Circle_lim_h02_P4_DPe0.380E+00_NR0002';
meshpath =  matmeshpath;


%% serial
pos = strfind(solution,'_P');
for i=1:10
    if strcmp(solution(pos+i),'_')
        pos = pos+i-1;
        break
    end
end

meshname = [solution(7:pos) '.h5'];
HDF5load([meshpath,meshname])
HDF5load([solpath,solution])
userial = transpose(reshape(u,[2,numel(u)/2])); 
utserial = transpose(reshape(u_tilde,[2,numel(u_tilde)/2])); 

if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end

refEl = createReferenceElement(elemType,size(T,2));

%% parallel
for ip = 1:nproc
    meshname = [solution(7:pos) '_' num2str(ip) '_' num2str(nproc) '.h5'];
    solname = [solution '_' num2str(ip) '_' num2str(nproc) '.h5'];
    HDF5load([meshpath,meshname])
    HDF5load([solpath,solname])
    uparall = transpose(reshape(u,[2,numel(u)/2]));
    utparall = transpose(reshape(u_tilde,[2,numel(u_tilde)/2]));
    
    ind = col(transpose(bsxfun(@plus,transpose((loc2glob_el-1)*size(refEl.NodesCoord,1)),(1:size(refEl.NodesCoord,1)))));
    indf = col(transpose(bsxfun(@plus,transpose((loc2glob_fa-1)*size(refEl.NodesCoord1d,1)),(1:size(refEl.NodesCoord1d,1)))));
    max(abs(userial(ind,:)-uparall))
    max(abs(utserial(indf,:)-utparall))

    
    
end