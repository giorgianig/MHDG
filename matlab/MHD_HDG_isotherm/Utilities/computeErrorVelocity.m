function [error eg] = computeErrorVelocity(M_star,u_star,u_int)

Ne = size(M_star,3);
Nv = size(M_star,1);
errVel = u_star-u_int;
errVel = transpose(reshape(errVel,2,Ne*Nv));
u_star = transpose(reshape(u_star,2,Ne*Nv));
error2 = zeros(Ne,1);
u2 = error2;
int_u02 = 0;
for iElem = 1:Ne
    
    ind = (iElem-1)*Nv + (1:Nv);
    Me = M_star(:,:,iElem);
    ex = errVel(ind,1);
    ey = errVel(ind,2);
    ux = u_star(ind,1);
    uy = u_star(ind,2);
%     error(ielem) = sqrt(ex'*Me*ex + ey'*Me*ey)/sqrt(ux'*Me*ux + uy'*Me*uy);
    
    error2(iElem) = (ex'*Me*ex + ey'*Me*ey);
    u2(iElem) = (ux'*Me*ux + uy'*Me*uy);
    int_u02 = int_u02 + (ux'*Me*ux + uy'*Me*uy);
end
    
error = sqrt(error2/int_u02);
eg = sqrt(sum(error2)/int_u02);



