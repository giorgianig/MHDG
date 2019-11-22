function L2err = computeErrorAnalyticSol(X,T,u,refEl)

global neq


% mesh informations
Ne = size(T,1);
Nv = size(T,2);

% reshape u
u = transpose(reshape(u,neq,Nv*Ne)); % Nv·Ne x 2

% allocate error vector
err2 = zeros(neq,1);
intsol = err2;

for iElem = 1:Ne
       
    % mesh
    Te = T(iElem,:);
    Xe = X(Te,:);
    
    % indices
    ind_u = (iElem-1)*Nv + (1:Nv);
    
    % elemental error and velocity
    ue = u(ind_u,:);
    
    % compute error
    [err2e, intsole] = computeElementalError(ue,Xe,refEl);
    err2 = err2 + err2e;
    intsol = intsol + intsole;
end
L2err = sqrt(err2./intsol);

function [err2, intsol] = computeElementalError(sol,Xe,refEl)

global neq

% reference element
IPw = refEl.IPweights;
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;

% number of Gauss points
ngauss = length(IPw);

% coordinates of gauss points
Xg = N*Xe;

% solution at gauss points
solg = N*sol;

% analytical solution at gauss points
exsol = analyticalSolution(Xg);

% error at gauss points
err_g = solg-exsol;

% initialization
err2 = zeros(neq,1);
intsol = err2;
for g = 1:ngauss
    
    % Values at current integration point
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    % Jacobian
    J = [Nxi_g*Xe(:,1)	  Nxi_g*Xe(:,2)   
         Neta_g*Xe(:,1)  Neta_g*Xe(:,2)];
     
    % Integration weight
    dvolu=IPw(g)*det(J);
    
    %Contribution of the current integration point to the elemental L2 Norm 
    err2 = err2 + transpose(err_g(g,:).^2)*dvolu;
    intsol = intsol + transpose(exsol(g,:).^2)*dvolu;
end



