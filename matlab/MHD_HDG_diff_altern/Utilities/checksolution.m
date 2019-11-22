function u = checksolution(u)

global neq useThreshold
u = transpose(reshape(u,neq,numel(u)/neq));
ind = u(:,1)<useThreshold;

if any(ind)
    disp('Using threshold')
end
u(ind,1) = useThreshold;  
u(ind,2) = 0;
u = col(u');
u = u(:);

