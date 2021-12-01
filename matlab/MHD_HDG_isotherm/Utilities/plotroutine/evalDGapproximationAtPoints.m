function up = evalDGapproximationAtPoints(xs,u,X,T,refEl)
% Ux=computeDGapproximationAtPoints(x,u,X,T,refEl)
% Evaluates the DG approximation given by the nodal values in u. u may have
% several columns (several field components)

X = X/1.901e-3;
xs = xs/1.901e-3;






np = size(T,2);
up = zeros(size(xs,1),size(u,2));
npoints = size(xs,1);

% tol = 1e-10;
tol = 1e-2;
correl = zeros(size(xs,1),1);
% Find corresponding elements
for iel = 1:size(T,1)
    Xe = X(T(iel,:),:);
    bc = [Xe(1:3,1)';Xe(1:3,2)';ones(1,3)]\[xs(:,1)';xs(:,2)';ones(1,npoints)];
    correl(all([bc>=-tol; bc<=1+tol],1)) = iel ;
    
end


% ratio = 0.1;
% %Loop in points
% for i=1:npoints
%     x = xs(i,:);
%     el = correl(i);
%     if ~el, continue,end
%     ind = (el-1)*np+(1:np);
% %     [e,xieta] = findElementForx(x,X,T,refEl);
% %     if(e==0) error('Error in evalDGapproximationAtPoints: x=(%g,%g) is not in the mesh',x(1),x(2)); end
% % %     if (e==0), continue, end
%     xieta = inverseIsopTransf(x,X(T(el,:),:),refEl);
%     shapeFunctions = computeShapeFunctionsAtPoints(refEl.degree,refEl.NodesCoord,xieta);
%     N = shapeFunctions(:,:,1)';
%     ue = u(ind,:);
%     up(i,:)= N*ue;
%     
%     if i>(npoints*ratio)
%         disp(['Done ', num2str(ratio*100) '%'])
%         ratio = ratio+0.1;
%     end
% end

%Loop in points
for el=unique(correl)'
    if ~el, continue,end
    loc = correl == el;
    x = xs(loc,:);
    ind = (el-1)*np+(1:np);
%     [e,xieta] = findElementForx(x,X,T,refEl);
%     if(e==0) error('Error in evalDGapproximationAtPoints: x=(%g,%g) is not in the mesh',x(1),x(2)); end
% %     if (e==0), continue, end
    xieta = inverseIsopTransf(x,X(T(el,:),:),refEl);
    shapeFunctions = computeShapeFunctionsAtPoints(refEl.degree,refEl.NodesCoord,xieta);
    N = shapeFunctions(:,:,1)';
    ue = u(ind,:);
    if any(ue(:,1)==0), stop,end
    up(loc,:)= N*ue;
%     if any(up(loc,:)<0),stop,end
    if any(isnan(N*ue)), stop, end
end
