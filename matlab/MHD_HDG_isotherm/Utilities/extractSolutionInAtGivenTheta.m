function upol = extractSolutionInAtGivenTheta(u,T,refElPol,refElTor,t)

global ntor theta

Np2d      = size(refElPol.NodesCoord,1);              % Number of points for 2d element
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nv     = Np1dTor*Np2d;                                      % Number of points for 3d element
N2d   = size(T,1);                                                % Number of 2d elements

% find toroidal element
te = linspace(0,theta,ntor+1);
if ismember(t,te)
    torel = find(te-t==0);
    torel=torel-1;
    if torel<1,torel=1;end
else
    d = te-t;
    torel = find(d(2:end).*d(1:end-1)<0);
end
tlim = [te(torel),te(torel+1)];
trel = -1+2*(t-tlim(1))/(tlim(2)-tlim(1));

% create 1d shape functions
V = Vandermonde_LP(refElTor.degree,refElTor.NodesCoord1d);
p = orthopoly1D_deriv(trel,refElTor.degree);
N = transpose(V'\p);

% interpolate at given theta
upol = zeros(N2d*Np2d,1);
for iel = 1:N2d
        iElem = (torel-1)*N2d+iel; % 3d numbering of the element
        ind_ue = (iElem-1)*Nv + (1:Nv);
        ind_upol = (iel-1)*Np2d + (1:Np2d);
        ue = u(ind_ue);
        ue = transpose(reshape(ue, Np2d,Np1dTor));        
        upol(ind_upol) = N*ue;
end

