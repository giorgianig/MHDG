function [physicalPoints]= getHOcurve (nDeg,listEdgeNodes,X,shapeFunctions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - nDeg = degree of the polynomial interpolation                         %
% - listEdgeNodes = set of nodes that define the edge of degree p that we %
%   want to have more points of                                           %
% - X = coordinates matrix                                                %
% - midPoints = [ empty ] -> the output of the function is a dense set of %
%                            points of the HO curve that defines de edge  %
%               [ 1 ] -----> the output of the function is the midPoints  %
%                            of the nodes of the edge computed in the     %
%                            reference element                            %
%               [ 2 ] -----> the output of the function is a new          %
%                            location of the edge nodes: an equispaced    %
%                            distribution on the reference element        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mapping of the reference points to the physical space

if(length(listEdgeNodes)>nDeg+1)
    %crec q no entra aqui mai... ho puc treure
    [physicalPoints]=computeImageHighOrderRepresentation(shapeFunctions,X(:,listEdgeNodes(1:2:(length(listEdgeNodes)))));
else   
    [physicalPoints]=computeImageHighOrderRepresentation(shapeFunctions,X(:,listEdgeNodes));
end

