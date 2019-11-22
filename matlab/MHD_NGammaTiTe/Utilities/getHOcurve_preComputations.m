function [shapeFunctions] = getHOcurve_preComputations(nDeg,midPoints,fekete_selection)

%% fekete distribution or equispaced
% if(nargin<5 || fekete_selection~=1)
%     fekete_selection=0;
% end 


%% Set the 1D reference element
ini_ref=-1;
fin_ref=1;
length_ref=fin_ref-ini_ref;

%% p+1 points that "define" the shape function

numBasePoints= nDeg+1;

if(fekete_selection==0)  
    %----  equispaced distribution  ----
    basePoints = (ini_ref  : length_ref/(numBasePoints-1) : fin_ref )';
else
    %----   feketee distribution   ----
    basePoints = feketeNodes1D(nDeg,1:numBasePoints);
end

% % la idea d'equispaiar a l'espai de referencia tal com al fisic no va be?...
% if(length(listEdgeNodes)>nDeg+1)
%     % aqui sembla q no hi entra res....
% %     aux = X(:,listEdgeNodes(1:2:(length(listEdgeNodes)-1)))-X(:,3:2:length(listEdgeNodes));
% else
%     aux = X(:,listEdgeNodes(1:nDeg))-X(:,listEdgeNodes((1:nDeg)+1));   
% end
% length_intervals = sqrt( sum(aux.*aux,1) );
% length_curve = sum(length_intervals,2);
% relative_length = length_intervals/length_curve*length_ref;
% 
% basePoints = zeros(numBasePoints,1);
% basePoints(1) = ini_ref;
% for i = 2:numBasePoints
%     basePoints(i) = basePoints(i-1) + relative_length(i-1);
% end



%% reference points in [-1;1] where we want the value of the shape function

if nargin<4 || midPoints==0 % used to plot polinomials as mini-lines
                            % output of nDeg*k points (a lot of points)
    k=10;
    numIntervals=nDeg*k;
    h=length_ref/numIntervals;
    referencePoints = ( ini_ref  : h : fin_ref )';
    
elseif midPoints==1 % used to bisect an edge: 
                    % output of 2*nDeg-1 points
    
    aux= (basePoints+fin_ref)/length_ref;
    referencePoints = [ aux+ini_ref  ;   aux(2:(nDeg+1))];

elseif midPoints==2 % used to equispace an edge
                    % output of p-1 points
    
    h=length_ref/(nDeg);
    referencePoints = ( (ini_ref+h)  : h : (fin_ref-h) )';

end

% Vandermonde matrix

V = Vandermonde_LP ( nDeg , basePoints );
[L,U,P] = lu(V');

%Shape function on the points over [-1,1]

p = orthopoly1D ( referencePoints  , nDeg ) ;   

N = U\(L\(P*p));

shapeFunctions = N;

