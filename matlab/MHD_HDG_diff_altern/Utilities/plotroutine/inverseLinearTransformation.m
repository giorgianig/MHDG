function xieta=inverseLinearTransformation(x,Xe)

x1=Xe(1,:)'; x2=Xe(2,:)'; x3=Xe(3,:)';


% if norm(x1-x')<1.e-10*norm(x')+1.e-14
%     xieta = [-1,1];
% elseif norm(x2-x')<1.e-10*norm(x')+1.e-14
%     xieta = [1,-1];
% elseif norm(x1-x')<1.e-10*norm(x')+1.e-14
%     xieta = [-1,1];
% else
    J = [ (x2-x1)/2 (x3-x1)/2 ];
%     xieta = ( J\( x'-(x2+x3)/2 ) )';
    auxx  = x(:,1)'-(x2(1)+x3(1))/2 ;
    auxy  = x(:,2)'-(x2(2)+x3(2))/2 ;
    xieta = ( J\( [auxx; auxy]) )';
%end


