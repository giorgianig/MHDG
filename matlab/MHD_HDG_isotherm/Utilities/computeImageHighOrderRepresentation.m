function [physicalPoints]=computeImageHighOrderRepresentation(shapeFunctions,physicalNodes)

%% Computation of the shape functions on the desired reference points

% % Vandermonde matrix
% 
% V = Vandermonde_LP ( nDeg , basePoints );
% [L,U,P] = lu(V');
% 
% %Shape function on the points over [-1,1]
% 
% p = orthopoly1D ( referencePoints  , nDeg ) ;   
% 
% N = U\(L\(P*p));
% 
% shapeFunctions = N;


%% Mapping of the reference points to the physical space

physicalPoints = physicalNodes*shapeFunctions ;