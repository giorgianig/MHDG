function varargout = plotMesh_HOElement_mod(X,T,linewidth,linecolor,linestyle)
%   option (optional): type 'plotNodes' to see the nodes' position on the
%                      ploted mesh, or type 'plotNodesNum' to see their global
%                      number.
%   nodesNum (necesary if option = 'plotNodesNum'): type 'all' to plot the
%                                                   global postion of all
%                                                   nodes, or enter a list
%                                                   array with the selected
%                                                   nodes.
np = size(T,2);
element.type = 'tri';
element.order = 0.5*(-3+sqrt(9+4*(2*np-2)));
element.distribution = 'fekete';%'equispaced'; %'fekete';
element.numNod = np;
faceNodes = faceNodes_aux(element);
[nOfFaces,nOfNodesPerFace] = size(faceNodes);
oFaceNodes = zeros(1,nOfFaces*(nOfNodesPerFace-1));
np = nOfNodesPerFace - 1;
aux = 1 - np;
aux2 = 0;
for iface = 1:nOfFaces
    aux = aux + np;
    aux2 = aux2 + np;
    oFaceNodes(aux:aux2) = faceNodes(iface,1:np);
end

%Conectivity for the faces
patchFaces = T(:,oFaceNodes);

%Plot mesh
% patchHandle = patch('Faces',patchFaces,'Vertices',X,'FaceColor','none','EdgeAlpha',1,'LineWidth',2);
% axis equal

% Shape funtions 1D
[p_order]=element.order;
feketeDistribution = strcmp(element.distribution,'fekete');
[shapeFunctions] = getHOcurve_preComputations(p_order,0,feketeDistribution);

hold on
for i=1:size(T,1)
    face1=patchFaces(i, (1:(np+1))            );
    [physicalPoints1]= getHOcurve (p_order,face1,X,shapeFunctions);
    face2=patchFaces(i, ((np+1):(2*np+1))     );
    [physicalPoints2]= getHOcurve (p_order,face2,X,shapeFunctions);
    if(nOfFaces == 3)
        face3=patchFaces(i,[  (2*np+1):(3*np) 1  ]      ) ;
        [physicalPoints3]= getHOcurve (p_order,face3,X,shapeFunctions);
        physicalPoints = [ physicalPoints1 physicalPoints2 physicalPoints3 ]';
    elseif(nOfFaces == 4)
        face3=patchFaces(i,[ (2*np+1):(3*np+1) ]      ) ;
        [physicalPoints3]= getHOcurve (p_order,face3,X,shapeFunctions);
        physicalPoints = [ physicalPoints1 physicalPoints2 physicalPoints3 ]';
        face4=patchFaces(i,[ (3*np+1):(4*np) 1 ]      ) ;
        [physicalPoints4]= getHOcurve (p_order,face4,X,shapeFunctions);
        physicalPoints = [physicalPoints; physicalPoints4'];
    end
    
    numPoints=size(physicalPoints,1);
    
    patchHandle = patch('Faces',1:numPoints,'Vertices',physicalPoints,...
        'FaceColor','none','EdgeAlpha',1,'LineWidth',linewidth,'Linestyle',linestyle,'edgecolor',linecolor);    
end

%% final configurations
axis equal
%     caxis([0 1])

%% Output variable
if ~nargout
    varargout = [];
else
       varargout = {patchHandle};
%     varargout = [];
    % no usem argments de sortida de moment
end
end

function res = faceNodes_aux(element)
Tref = 1:element.numNod;
switch element.type
    case 'tri'
        res = [...
            getEdge_tri(Tref,element.order,1)
            getEdge_tri(Tref,element.order,2)
            getEdge_tri(Tref,element.order,3)
            ];
    case 'quad'
        res = [...
            getEdge_quad(Tref,element.order,1)
            getEdge_quad(Tref,element.order,2)
            getEdge_quad(Tref,element.order,3)
            getEdge_quad(Tref,element.order,4)
            ];
    otherwise
        error('Cannot plot this type of element yet')
end
end