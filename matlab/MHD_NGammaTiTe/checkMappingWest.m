matmeshpath = '/home/giorgio/Matlab/Meshes/Meshes_2D/WEST/'; 
matmesh = load([matmeshpath,meshname(1:end-3) '.mat']);
F = hdg_Preprocess(matmesh.T,matmesh.elementFaceInfo,{'Tb_IN','Tb_OUT'},refEl,[0 0]);

iel = 31664;
Fe = F(iel,:);
np = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
neq = 4;

ind_ue = col((iel-1)*np*neq + (1:np*neq));
ind_ut = col(bsxfun(@plus,(Fe-1)*nv*neq,(1:nv*neq)'));


ufort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/solu.txt');
utildefort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/solutilde.txt');
UU = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/UU.txt');
U0 = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/U0.txt');

ue = ufort(ind_ue);
uf  = utildefort(ind_ut);
norm(UU*uf+U0-ue)

ue = transpose(reshape(ue,neq,numel(ue)/neq));
uf  = transpose(reshape(uf,neq,numel(uf)/neq));