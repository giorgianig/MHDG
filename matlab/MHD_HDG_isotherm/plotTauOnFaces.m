clear
% close all
lscale = 1.901*1e-3;

ms = 20;

meshpath = '/home/giorgio/Dropbox/Fortran/MHDG/test/test_Temperature/West/';
meshname = 'West_Hugo_h0.02_refBound0.01_QUADS_P3.h5';
tau = readMatTxt([meshpath,'tau_save.txt']);
xy = readMatTxt([meshpath,'xy_g_save.txt']);
tauf = readMatTxt([meshpath,'tau_save_bound.txt']);
xyf = readMatTxt([meshpath,'xy_g_save_bound.txt']);
elemType = 0;


tau = [tau;tauf];
xy = [xy;xyf];

% tau = tauf;
% xy = xyf;

clear tauf xyf

indzer = sum(xy,2)==0;
tau(indzer,:) = [];
xy(indzer,:) = [];


HDF5load([meshpath,meshname])
if elemType==1
    elemType=0;
elseif elemType==0
    elemType=1;
end
% X(:,1) = X(:,1) +3.4;

X = X/lscale;


x = xy(:,1);
y = xy(:,2);

%% plot tau 1
figure
plotMesh(X,T,elemType)
C = tau(:,1);

cdivs = 100;
[~, edges] = hist(C,cdivs-1);
edges = [-Inf edges Inf]; % to include all points
[Nk, bink] = histc(C,edges);

hold on;
cmap = jet(cdivs);
for ii=1:cdivs
    idx = bink==ii;
    plot(x(idx),y(idx),'.','MarkerSize',ms,'Color',cmap(ii,:));
end

colormap(cmap)
if (min(C)~=max(C))
    caxis([min(C) max(C)])
else 
    caxis([min(C)-1 max(C)+1])
end
colorbar
title('tau on 1st equation')





%% plot tau 2
figure
plotMesh(X,T,elemType)
C = tau(:,2);

cdivs = 100;
[~, edges] = hist(C,cdivs-1);
edges = [-Inf edges Inf]; % to include all points
[Nk, bink] = histc(C,edges);

hold on;
cmap = jet(cdivs);
for ii=1:cdivs
    idx = bink==ii;
    plot(x(idx),y(idx),'.','MarkerSize',ms,'Color',cmap(ii,:));
end

colormap(cmap)
if (min(C)~=max(C))
    caxis([min(C) max(C)])
else 
    caxis([min(C)-1 max(C)+1])
end
colorbar
title('tau on 2nd equation')

