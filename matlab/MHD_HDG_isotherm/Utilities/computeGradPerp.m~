function res  = computeGradPerp(X)

global Mesh

[~,uxex,uyex] = analyticalSolution(X);
b = defineMagneticField(X);
GradU = ([permute(uxex,[2 3 1]),permute(uyex,[2 3 1])]);
GradPerpU = zeros(size(GradU));
for ip = 1:size(Mesh.X,1)
    bGradU = GradU(:,:,ip)*transpose(b(ip,:));
    GradPerpU(:,:,ip) = GradU(:,:,ip)-bGradU*b(ip,:);
end

res = col(permute(GradPerpU,[2 1 3]));
