function intsol = computeIntSol(X,T,u,refEl)

global neq

if size(u,2)==1
    stop
end
% mesh informations
Ne = size(T,1);
Nv = size(T,2);

% reshape u
% u = transpose(reshape(u,neq,Nv*Ne)); % Nv·Ne x neq

% allocate error vector
intsol = zeros(Ne,neq);

for iElem = 1:Ne
       
    % mesh
    Te = T(iElem,:);
    Xe = X(Te,:);
    
    % indices
    ind_u = (iElem-1)*Nv + (1:Nv);
    
    % elemental error and velocity
    ue = u(ind_u,:);
    
    % compute error
    intsole = computeElementalError(ue,Xe,refEl);
    intsol(iElem,:) = intsole;
end

function intsol = computeElementalError(sol,Xe,refEl)

global neq

% reference element
IPw = refEl.IPweights;
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;

% number of Gauss points
ngauss = length(IPw);

% solution at gauss points
solg = N*sol;

xg = N*Xe(:,1);
% initialization
intsol = zeros(neq,1);
for g = 1:ngauss
    
    % Values at current integration point
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    % Jacobian
    J = [Nxi_g*Xe(:,1)	  Nxi_g*Xe(:,2)   
         Neta_g*Xe(:,1)  Neta_g*Xe(:,2)];
     
    % Integration weight
    dvolu=IPw(g)*det(J)*xg(g);
    
    %Contribution of the current integration point to the elemental L2 Norm 
    intsol(:) = intsol(:) + transpose(solg(g,:))*dvolu;
end



