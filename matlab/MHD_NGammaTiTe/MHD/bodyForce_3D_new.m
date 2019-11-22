function f = bodyForce_3D_new(X,t)

global testcase

switch testcase.n
    case 1
        f = bf1_3D(X,t);
    case 2
        f = bf2_3D(X,t);
end

