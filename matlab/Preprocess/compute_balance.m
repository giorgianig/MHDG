function [bal,nel] = compute_balance(x,xb)
% 
% xex = [0,x,2*pi];
% ndiv = numel(x)+1;
% nel = zeros(ndiv,1);
% for it=1:ndiv
% %     ind = all([xb>xex(it) xb<xex(it+1)],2);
% %       ind=bsxfun(@minus,xb,xex(1:end-1)).*bsxfun(@minus,xb,xex(2:end))<0;
%     ind = (xb-xex(it)).*(xb-xex(it+1))<0;
%     nel(it) = nnz(ind);
% end
% % nel=sum(ind,1);
% nel_mean=sum(nel)/ndiv;
% % f = sum((nel-nel_mean).^2);
ndiv = numel(x)-1;
nel = zeros(ndiv,1);
for it=1:ndiv
    ind = all([xb>x(it) xb<x(it+1)],2);
    nel(it) = nnz(ind);
end
nel_mean=sum(nel)/ndiv;
bal = (nel-nel_mean)/nel_mean;