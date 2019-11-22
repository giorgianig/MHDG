function varargout = plotMesh_HOElement(X,T,option,input)%nodesNum,qualityElements)
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

        %% colouring    
        if(nargin>2)
            for iOption = 1:length(option)
                hold on 
                switch option{iOption}
                    case 'plotHOquality'
                        patch(  physicalPoints(1:(size(physicalPoints,1)),1) ,...
                                physicalPoints(1:(size(physicalPoints,1)),2) ,...
                                input.quality(i),...
                                'EdgeAlpha',1 ) %'LineWidth',1,...
                    case 'plotHOquality3D'
                        patch(  physicalPoints(1:(size(physicalPoints,1)),1) ,...
                                physicalPoints(1:(size(physicalPoints,1)),2) ,...
                                physicalPoints(1:(size(physicalPoints,1)),3) ,...
                                input.quality(i),...
                                'EdgeAlpha',1 ) %'LineWidth',1,...
                    case 'plotSurfaces'
                        patch(  physicalPoints(1:(size(physicalPoints,1)),1) ,...
                                physicalPoints(1:(size(physicalPoints,1)),2) ,...
                                physicalPoints(1:(size(physicalPoints,1)),3) ,...
                                input.quality,...
                                'EdgeAlpha',1 ) %'LineWidth',1,...
                    case 'fatLines'
                        patch('Faces',1:numPoints,'Vertices',physicalPoints,'FaceColor','none','EdgeAlpha',1,'LineWidth',2);        
                    case 'green'
                        patch('Faces',1:numPoints,'Vertices',physicalPoints,'FaceColor','none','EdgeAlpha',1,'LineWidth',1,'EdgeColor',[0 0.7 0]);        
                    case {'red', 'blue', 'black'}
                        patch('Faces',1:numPoints,'Vertices',physicalPoints,'FaceColor','none','EdgeAlpha',1,'LineWidth',1,'EdgeColor',option{iOption});  
                    otherwise
                        patch('Faces',1:numPoints,'Vertices',physicalPoints,'FaceColor','none','EdgeAlpha',1,'LineWidth',1);
                end
            end
        else
            patch('Faces',1:numPoints,'Vertices',physicalPoints,...
                'FaceColor','none','EdgeAlpha',1,'LineWidth',1,'Linestyle','-','edgecolor','k');
        end

        %%
    %     hold on
    %     plot(physicalPoints1(1,:),physicalPoints1(2,:),'o');
    %     plot(physicalPoints2(1,:),physicalPoints2(2,:),'r*');
    %     plot(physicalPoints3(1,:),physicalPoints3(2,:),'g+');
    %     %%
    end
% plot(X(:,1),X(:,2),'ko','markersize',15);
% axis equal

%    list = 1:size(X,1);
%    fontSize = 10;
% for inode = list
%    text(X(inode,1),X(inode,2),X(inode,3),int2str(inode),'FontSize',fontSize,...
%        'Color',[1 0 0])
% end

    %% Optional plots: num nodes, elements...
    if(nargin>2)
        X = X';
        for iOption = 1:length(option)
            hold on 
            switch option{iOption}
                case 'plotNodes'
%                     plot(X(:,1),X(:,2),'o','markerSize',5,'markerFaceColor','n','markerEdgeColor','black')
                        plot(X(:,1),X(:,2),'*','markerSize',10,'markerFaceColor','b','markerEdgeColor','b')

                case 'plotNodes3D'
%                     plot3(X(:,1),X(:,2),X(:,3),'o','markerSize',12,'markerFaceColor','n','markerEdgeColor','black')
                        plot3(X(:,1),X(:,2),X(:,3),'*','markerSize',12,'markerFaceColor','n','markerEdgeColor','black')

                case 'plotNodesNum'
                   if(isfield(input,'nodesNum') == false)
                       list = 1:size(X,1);
                       fontSize = 10;
                   else
                       list = input.nodesNum;
                       fontSize = 15;
                       plot(X(list,1),X(list,2),'o','markerSize',3,'markerFaceColor','b')
                   end
                   for inode = list
                       text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize,...
                           'Color',[1 0 0])
                   end

                case 'plotNodesNumAndElements'
                   if(nargin < 5 || isfield(input,'nodesNum') == false)
                       list = 1:size(X,1);
                       fontSize = 10;
                   else
                       list = input.nodesNum;
                       fontSize = 15;
                       plot(X(list,1),X(list,2),'o','markerSize',3,'markerFaceColor','b')
                   end
                   for inode = list
                       text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize,...
                           'Color',[1 0 0])
                   end
                   for iElem = 1:size(T,1)
                       xbar = 1/3*(X(T(iElem,1),1)+X(T(iElem,2),1)+X(T(iElem,3),1));
                       ybar = 1/3*(X(T(iElem,1),2)+X(T(iElem,2),2)+X(T(iElem,3),2));
                       text(xbar,ybar,int2str(iElem),'FontSize',fontSize+2,...
                           'Color',[0 0 1])
                   end

                case 'plotElements'
                   fontSize =16;
                   for iElem = 1:size(T,1)
                       xbar = 1/3*(X(T(iElem,1),1)+X(T(iElem,2),1)+X(T(iElem,3),1));
                       ybar = 1/3*(X(T(iElem,1),2)+X(T(iElem,2),2)+X(T(iElem,3),2));
                       text(xbar,ybar,int2str(iElem),'FontSize',fontSize,...
                           'Color',[0 0 1])
                   end
            end
        end
       hold off
    end
    
    %% final configurations
    axis equal
%     caxis([0 1])

    %% Output variable
    if ~nargout
       varargout = [];
    else
    %    varargout = {patchHandle};
    varargout = [];
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