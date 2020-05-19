clear

p =2% [1,5,8,10];
meshName = 'Circle_CAMILLE';
path2Mesh = '/home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/';

for ip = p
    
    load([path2Mesh meshName '_P' num2str(ip) '.mat'])
    aux = whos('Tb_*');
    boundaryNames = {aux.name};
    fixDataWithRepeatedNodes
    save([path2Mesh meshName '_P' num2str(ip) '_fix.mat'],'X','T',...
        'elementFaceInfo','elemInfo','Tb_*')
    clear Tb_* X T elemInfo elemFaceInfo

    
end 