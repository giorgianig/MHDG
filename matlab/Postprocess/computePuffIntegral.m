function [puff,area] = computePuffIntegral(X,T,q,refEl,extfaces,boundaryFlags,simpar)

neq = simpar.Neq;
ndm = simpar.Ndim;
qres = transpose(reshape(q,[neq*ndm,numel(q)/ndm/neq]));
ind_puff = boundaryFlags==6;

nf = size(extfaces,1);
np = size(refEl.N,2);
puff =0;
area = 0;
for i=1:nf
    
    if (~ind_puff(i)),continue,end
    
    el = double(extfaces(i,1));
    fa = double(extfaces(i,2));    
    ind=(el-1)*np + (1:np);
    indf = ind(refEl.faceNodes(fa,:));
    Xf = X(T(el,refEl.faceNodes(fa,:)),:);
    qf = qres(indf,(neq-1)*ndm+1:neq*ndm);
    [inte,areae] = computeElementalIntegral(Xf,qf,refEl);
    puff = puff+inte;
    area = area+areae;
    
end

puff = -puff*simpar.physics.diff_nn*simpar.adimensionalization.diffusion_scale*...
    simpar.adimensionalization.density_scale/simpar.adimensionalization.length_scale;

function [inte,area]=computeElementalIntegral(Xf,qf,refEl)


ngauss = size(refEl.IPweights1d);
N = refEl.N1d;
Nx = refEl.N1dxi;
xyg = N*Xf;
xyg_p = Nx*Xf;
qfg = N*qf;
inte = 0;
area = 0;
for g=1:ngauss
        
    % Integration weight
    xyDerNorm_g = norm(xyg_p(g,:));
    dline = refEl.IPweights1d(g)*xyDerNorm_g;
    dline = 2*pi*dline*xyg(g,1);
        
    % Unit normal to the boundary
    t_g = xyg_p(g,:)/xyDerNorm_g;
    n_g = [t_g(2) -t_g(1)];
    
    qdotn = dot(qfg(g,:),n_g);
    inte = inte + qdotn*dline;
    area = area+dline;

end