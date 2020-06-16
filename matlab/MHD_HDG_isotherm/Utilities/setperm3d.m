function res = setperm3d(N1dPol,N1dTor,Neq)
% check
res = colt(     flipud(bsxfun(@plus,(1:N1dPol)',( (1:N1dTor)-1)*N1dPol   ) )     );
res = col( bsxfun(@plus,(1:Neq)',(res-1)*Neq )); 