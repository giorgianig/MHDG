function [e]=getEdge_tri(T,p,localEdge)

if(localEdge==1)
    [e] = T( : , [1  4:(3+(p-1))              2] ) ;
elseif(localEdge==2)
    [e] = T( : , [2  (4+(p-1)):(3+2*(p-1))    3] ) ;
else
    [e] = T( : , [3  (4+2*(p-1)):(3+3*(p-1))  1] ) ;
end

% if(localEdge==1)
%     [e] = getFirstEdge(Telem,p) ;
% elseif(localEdge==2)
%     [e] = getSecondEdge(Telem,p) ;
% else
%     [e] = getThirdEdge(Telem,p) ;
% end