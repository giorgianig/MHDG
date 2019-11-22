clear

p =3% [1,5,8,10];
meshName = 'West_Hugo_h0.04_refCorn';
path2Mesh = '/home/giorgio/Matlab/Meshes/Meshes_2D/WEST/';

for ip = p
    
    load([path2Mesh meshName '_P' num2str(ip) '.mat'])
    aux = whos('Tb_*');
    boundaryNames = {aux.name};
    fixDataWithRepeatedNodes
%     save([path2Mesh meshName '_P' num2str(ip) '_fix.mat'],'X','T',...
%         'elementFaceInfo','elemInfo','Tb_*')
%     clear Tb_* X T elemInfo elemFaceInfo

    
end 