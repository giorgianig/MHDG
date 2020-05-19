el = 48;

plotMesh(X,T)
hold on

for i = 1:size(T,2)
    text(X(T(el,i),1),X(T(el,i),2),num2str(T(el,i)),'color','g');
end