function u = checksolution(u,X,T,refEl)

global neq useThreshold physValueGar Mref
np = size(T,2);
ne = size(T,1);

if ~physValueGar
    u = transpose(reshape(u,neq,numel(u)/neq));
    
    ind_rho = u(:,1)<useThreshold;
    ind_eni = u(:,3)<useThreshold;
    ind_ene = u(:,4)<useThreshold;
    
    if any(ind_rho) || any(ind_eni) || any(ind_ene)
        disp('Using threshold')
    end
    u(ind_rho,1) = useThreshold;
    u(ind_rho,2) = 0;
    u(ind_eni,3) = useThreshold;
    u(ind_ene,4) = useThreshold;
    u = col(u');
    u = u(:);
else
    u = transpose(reshape(u,neq,numel(u)/neq));
    u1 = u(:,1);
    u2 = u(:,2);
    u3 = u(:,3);
    u4 = u(:,4);    
    U1 = col(refEl.N*reshape(u1,[np,ne])); % ng * nel: U1 at Gauss points in the elements
    U2 = col(refEl.N*reshape(u2,[np,ne])); % ng * nel: U2 at Gauss points in the elements
    U3 = col(refEl.N*reshape(u3,[np,ne])); % ng * nel: U3 at Gauss points in the elements
    U4 = col(refEl.N*reshape(u4,[np,ne])); % ng * nel: U4 at Gauss points in the elements
    up = cons2phys([U1,U2,U3,U4]);
%     pi = reshape(up(:,5),[ng,nel]);
%     pe = reshape(up(:,6),[ng,nel]);
%     ind_i = any(pi<useThreshold,1);
%     ind_e = any(pe<useThreshold,1);
   pi = up(:,5);
   pe = up(:,6);
   minpi = min(pi);
   minpe = min(pe);
   
   if minpi<useThreshold
       u3  = u3+1.5*Mref*abs(minpi);
       disp('Increasing ions energy')
   end
   if minpi<useThreshold
       u4  = u4+1.5*Mref*abs(minpe);
       disp('Increasing electrons energy')
   end
  u = col([u1,u2,u3,u4]');
end

