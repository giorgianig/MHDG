function [elementNum,xieta] = findElementForx(x,X,T,refEl)

tol = 1.e-3;


distance = sqrt((x(1)-X(:,1)).^2+(x(2)-X(:,2)).^2);
[mind,i]=min(distance); %node with minimum distance to x

[elements,j]=find(T==i);

elementNum = 0;

for i=1:length(elements)
    e = elements(i);
    Xe = X(T(e,:),:);
    cu = convhull(Xe(:,1),Xe(:,2));
    inout = inpolygon(x(1),x(2),Xe(cu,1),Xe(cu,2));
%     bc = [Xe(1:3,1)';Xe(1:3,2)';ones(1,3)]\[x(1);x(2);1];
    if ~inout, continue,end
%     xieta = inverseIsoparametricTransformation(x,Xe,refEl.NodesCoord,refEl.degree);
    xieta = inverseIsopTransf(x,Xe,refEl);
    if all(xieta<=1+tol) &  all(xieta>=-1-tol) & (xieta(1)+xieta(2))<=tol
       elementNum = e; 
       break;
    end
end

if elementNum == 0 %the point is not in the star of x_i -> search at neighbour elements...
    d = radiousStarNode(X,T(elements,:));
    iCloseNodes = find(distance<=d);
    elements2 = [];
    for i=1:length(iCloseNodes)
        [e,j]=find(T==iCloseNodes(i));
        elements2 = union(elements2,e);
    end
    elements = setdiff(elements2,elements);
    for i=1:length(elements)
        e = elements(i);
        Xe = X(T(e,:),:);
%         xieta = inverseIsoparametricTransformation(x,Xe,refEl.NodesCoord,refEl.degree);
        xieta = inverseIsopTransf(x,Xe,refEl);
        if all(xieta<=1+tol) &  all(xieta>=-1-tol) & (xieta(1)+xieta(2))<=tol
            elementNum = e;
            break;
        end
    end
    
end


%__________________________________________________________________
function r=radiousStarNode(X,T)
r = 0;
for e=1:size(T,1)
    Xe = X(T(e,:),:);
    H = max([norm(Xe(1,:)-Xe(2,:)),norm(Xe(2,:)-Xe(3,:)),norm(Xe(3,:)-Xe(1,:))]); 
    r = max(H,r);
end

