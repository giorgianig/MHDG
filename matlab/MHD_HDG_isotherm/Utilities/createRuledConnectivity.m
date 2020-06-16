function ruledConnectivity = createRuledConnectivity(X,T)

rule = 'x.^2 + y.^2 <= 2'; % T shape
disp(['ruled connectivity used:' rule])

X1 = X(:,1);
X2 = X(:,2);

x = X1(T);
y = X2(T);

eval(['auxmatrix = ' rule ';'])

ruledConnectivity = all(auxmatrix,2);