function [cf xyprof] = thwaites(X,T,u,refEl,infoFaces,Re)

nx = 200;
ny = 100;

[delta Ue xybl xyprof] = boundaryLayerProperties(X,T,u,refEl,infoFaces,3,nx,ny,100);

dx = sqrt((xyprof(2:end,1)-xyprof(1:end-1,1)).^2 + (xyprof(2:end,2)-xyprof(1:end-1,2)).^2);
ind = dx<1e-6;
dx = dx(~ind);
xyprof = xyprof(~ind,:);
Ue = Ue(~ind);

% smooth Ue
slot = 2*ceil(numel(Ue)/40);
for i= 1:numel(Ue)
    if i<slot/2
        delta = 1:2*i-1;
    elseif i>(numel(Ue)-slot/2)
        aux = numel(Ue)-i;
        delta = i-aux:numel(Ue);
    else
        delta = i+1-slot/2:i+slot/2;
    end
    Ue(i) = mean(Ue(delta));
end
    
d = cumsum(dx)-dx(1);

% dUe/dx
dudx = diff(Ue)./diff(d);
dudx = [dudx;0];

% integral Ue
intUe = zeros(size(d));
for i = 2:numel(d)
    intUe(i) = trapz(d(1:i),Ue(1:i).^5);
end

% theta
theta = sqrt(0.45/Re*intUe./Ue.^6);

% lambda
lambda = theta.^2*Re.*dudx;
indS = lambda>0.25;
indI = lambda<-0.09;
ind = indS+indI;

% H and L
H = fH(lambda);
L = fL(lambda);

% cf
cf = 2*L./(Re*theta);

cf = cf(~ind);
xyprof = xyprof(~ind);

function H = fH(lambda)
H = zeros(size(lambda));
H(lambda<0) = 2.088 + 0.0731./(lambda(lambda<0)+0.14);
H(lambda >= 0) = 2.61 - 3.75*lambda(lambda >= 0) + 5.24*lambda(lambda >= 0).^2;

function L = fL(lambda)
L = zeros(size(lambda));
L(lambda<0) = 0.22 + 1.402*lambda(lambda<0) +(0.018*lambda(lambda<0))./(lambda(lambda<0)+0.107);
L(lambda >= 0) = 0.22 + 1.57*lambda(lambda >= 0) - 1.8*lambda(lambda >= 0).^2;







