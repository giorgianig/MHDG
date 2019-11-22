function f = bodyForce_3D(X,t)

global testcase neq

N2d = size(X,1);
N1d = size(t,1);

% allocation
f = zeros(N2d*N1d,neq);

switch testcase.n
    case 1
        f = bf1_3D(X,t);
    case 2
        f = bf2_3D(X,t);
end

