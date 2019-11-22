function res = computeSourceIntgral(X,T,refEl)

% mesh informations
Ne = size(T,1);

% allocate error vector
res = 0;

for iElem = 1:Ne
    
    % mesh
    Te = T(iElem,:);
    Xe = X(Te,:);

    % compute integral
    integ_el = computeElementalError(Xe,refEl,iElem,Te);
    res = res + integ_el;
    
end

function res = computeElementalError(Xe,refEl,iElem,Te)

global testcase Magnetic

% reference element
IPw = refEl.IPweights;
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;

% number of Gauss points
ngauss = length(IPw);

% coordinates of gauss points
Xg = N*Xe;

% initialization
res = 0;
for g = 1:ngauss
    
    % Values at current integration point
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    % Jacobian
    J = [Nxi_g*Xe(:,1)	  Nxi_g*Xe(:,2)
        Neta_g*Xe(:,1)  Neta_g*Xe(:,2)];
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    
    % Source in the gauss point
    f = 0;
    switch testcase.n
        
        case 51 % Big hole
            if (Magnetic.flux2D(g,iElem)<=-0.88 && Magnetic.flux2D(g,iElem)>=-0.90)
                f = 3.20119388718018e-05;
            end
        case 52 % Small hole
            if (Magnetic.flux2D(g,iElem)<=-0.90 && Magnetic.flux2D(g,iElem)>=-1)
                f = 9.45155008295538e-06;
            end
        case 53 % No hole 
            if Magnetic.flux2D(g,iElem)<=-0.90
                f = 7.24032211339971e-06;
            end
        case 54 % No hole 
            if (Magnetic.flux2D(g,iElem)<=-0.98  && Magnetic.flux2D(g,iElem)>=-1)
                f = 5.78723047118297e-05;
            end        
        case 55 % No hole 
            if Magnetic.flux2D(g,iElem)<=-1.03
                f = 0.000115575293741846;
            end     
        case 56 % No hole - moving equilibrium
            phig = dot(N(g,:),Magnetic.flux2D(Te));
            if phig<=1.3 && phig<=1.2
                f = 3.811677788989013e-06;
            end                 
        otherwise
            f = bodyForce(Xg(g,:)); 
    end
    
    res = res + f*dvolu;
end



