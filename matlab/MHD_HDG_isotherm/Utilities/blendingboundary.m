function W = blendingboundary(Z,m,porder,s)
%
% W = blendingboundary(Z,m,porder,s)
%

if m == 1
    % First vertex ([0,0] in xi-eta)
    W = 1-Z(:,1)-Z(:,2);
elseif m == porder+1
    % Second vertex ([1,0] in xi-eta)
    W = Z(:,1);
elseif m == 2*porder+1
    % Third vertex ([0,0] in xi-eta)
    W = Z(:,2);
else
    % Edge nodes
    C = ones(porder-1,porder-1);
    C(:,1) = s(2:end-1).*(1-s(2:end-1));
    for i=2:porder-1
        C(:,i) = C(:,i-1).*s(2:end-1); 
    end
    C = inv(C);
    if m < porder+1 
        % First edge 
        ZC = Z; 
        ind = m-1;
    elseif m < 2*porder+1 
        % Second edge 
        ZC = [Z(:,2),1-Z(:,1)-Z(:,2)]; 
        ind = m-porder-1;
    else 
        % Third edge 
        ZC = [1-Z(:,1)-Z(:,2),Z(:,1)]; 
        ind = m-2*porder-1; 
    end
    W = Z(:,1)*0;
    for i=1:porder-1
        W = W + C(i,ind)*ZC(:,1).^i; 
    end
    W = W.*(1-ZC(:,1)-ZC(:,2));
end