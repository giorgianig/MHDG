function [shapeFun2d,gw2d,gp2d] = createShapeFunctions2dTensor(shapeFun1d,gw1d,gp1d,perm)

n1d = size(shapeFun1d,1);
ng  = size(shapeFun1d,2);

ind_x = reshape(gp1d,[ng,1]);
ind_y = reshape(gp1d,[1,ng]);
% aux_x = ind_x(:,ones(size(gp1d)),ones(size(gp1d)));
% aux_y = ind_y(ones(size(gp1d)),:,ones(size(gp1d)));
% gp2d  = [aux_x(:),aux_y(:)];

gp2d = [col(repmat(ind_x,1,numel(gp1d))) col(repmat(ind_y,numel(gp1d),1)) ];

gw1d_i = reshape(gw1d,[ng 1]);
gw1d_j = reshape(gw1d,[1 ng]);

gw2d = bsxfun(@times,gw1d_i,gw1d_j);
gw2d = gw2d(:);

shapeFun1d_i = reshape(shapeFun1d,[n1d 1 ng 1 2]);
shapeFun1d_j = reshape(shapeFun1d,[1 n1d 1 ng 2]);
shapeFun2d_1 = reshape(bsxfun(@times,shapeFun1d_i(:,:,:,:,1),shapeFun1d_j(:,:,:,:,1)),[n1d n1d ng^2]);
shapeFun2d_2 = reshape(bsxfun(@times,shapeFun1d_i(:,:,:,:,2),...
                                     shapeFun1d_j(:,:,:,:,1)),[n1d n1d ng^2]);
shapeFun2d_3 = reshape(bsxfun(@times,shapeFun1d_j(:,:,:,:,2),...
                                     shapeFun1d_i(:,:,:,:,1)),[n1d n1d ng^2]);
shapeFun2d_vec = cat(3,shapeFun2d_1,shapeFun2d_2,shapeFun2d_3);
shapeFun2d     = reshape(shapeFun2d_vec,[n1d^2,ng^2,3]);

% shapeFun2d = shapeFun2d([1,n1d,n1d^2,(n1d-1)*n1d+1,2:n1d-1,n1d+1:(n1d-1)*n1d,(n1d-1)*n1d+2:n1d^2-1],:,:);
% shapeFun2d = shapeFun2d(perm,:,:);
shapeFun2d(perm,:,:) = shapeFun2d;

% aux_1d = -1:2/(n1d-1):1;
% ind_x = reshape(aux_1d,[n1d,1,1]);
% ind_y = reshape(aux_1d,[1,n1d,1]);
% aux_x = ind_x(:,ones(size(aux_1d)));
% aux_y = ind_y(ones(size(aux_1d)),:);
% aux_2d  = [aux_x(:),aux_y(:)];
% 
% plot(aux_2d(:,1),aux_2d(:,2),'ro','markersize',10)
% for inode = 1:size(aux_2d,1)
%     text(aux_2d(inode,1),aux_2d(inode,2),num2str(inode),'fontsize',16)
% end
% 
% % aux_2d = aux_2d([1,n1d,n1d^2,(n1d-1)*n1d+1,2:n1d-1,n1d+1:(n1d-1)*n1d,(n1d-1)*n1d+2:n1d^2-1],:);
% % aux_2d = aux_2d(perm,:);
% aux_2d(perm,:) = aux_2d;
% 
% figure
% plot(aux_2d(:,1),aux_2d(:,2),'ro','markersize',10)
% for inode = 1:size(aux_2d,1)
%     text(aux_2d(inode,1),aux_2d(inode,2),num2str(inode),'fontsize',16)
% end
% stop