function res = filterSolution(T,u,refEl)

res = u;
thresh = 1e-3;
nel = size(T,1);
neq = 4;
coordRef = refEl.NodesCoord;
np = size(coordRef,1);
nDeg = refEl.degree;
npm1 = nDeg*(nDeg+1)/2;
V = Vandermonde_LP(nDeg,coordRef);
invV = V\eye(np);
s0 = log10(1/nDeg^4);

for iel = 1:nel
    
    ind = (iel-1)*np+(1:np);
    
    for ieq = 1:neq
        um = invV*u(ind,ieq);
        if sum(um.^2,1)>thresh
            continue
        end
        umho = zeros(size(um));
        umho(npm1+1:np,:) = um(npm1+1:np,:);
        se = log10( sum(umho.^2,1)./sum(um.^2,1))';
        if se>s0
            ulo = zeros(size(um));
            ulo(1) = um(1);
            res(ind,ieq) = V*ulo;
        end
    end
end
