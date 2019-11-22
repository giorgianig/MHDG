function h = findElementSize(X,T)

% finds the smallest face in all the elements

node1 = [X(T(:,1),1), X(T(:,1),2)]; % node 1 for all the elements in T
node2 = [X(T(:,2),1), X(T(:,2),2)]; % node 2 for all the elements in T
node3 = [X(T(:,3),1), X(T(:,3),2)]; % node 3 for all the elements in T

edge1 = node1-node2; % coord of first face (in the origin)
edge2 = node2-node3; % coord of second face (in the origin)
edge3 = node1-node3; % coord of third face (in the origin)

l1 = sqrt(edge1(:,1).^2+edge1(:,2).^2); % lenght of face 1 of all el in T
l2 = sqrt(edge2(:,1).^2+edge2(:,2).^2); % lenght of face 2 of all el in T
l3 = sqrt(edge3(:,1).^2+edge3(:,2).^2); % lenght of face 3 of all el in T

h = min([l1,l2,l3],[],2);
