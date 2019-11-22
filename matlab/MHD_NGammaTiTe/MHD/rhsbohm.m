function f = rhsbohm(X)

global testcase Mesh diff_pari diff_pare neq  epn Gmbohmi Gmbohme
global Mref
kpari = diff_pari;
kpare = diff_pare;
gmi = Gmbohmi;
gme = Gmbohme;



% allocation
f = zeros(size(X,1),neq);

% assign
x = X(:,1);
y = X(:,2);



switch testcase.n
    case 5
        a = pi/2;
        for i = 1:numel(x)
            xx = x(i);
            f(:,3) = ((0.11e1 - sin(a * xx)) * (0.2e1 - sin(a * xx)) + 0.2e1 / 0.3e1 * (1 - gmi) * (0.11e1 - sin(a * xx)) * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1)) * sin(a * xx) - kpari * (0.2e1 / 0.3e1 * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref) ^ epn * (-cos(a * xx) * a - cos(a * xx) * a * sin(a * xx)) / Mref / 0.45e2 - (0.11e1 - sin(a * xx)) * sin(a * xx) ^ 3 / 0.2e1;
            f(:,4) = ((0.11e1 - sin(a * xx)) * (0.14e1 - sin(a * xx)) + 0.2e1 / 0.3e1 * (1 - gme) * (0.11e1 - sin(a * xx)) * (0.14e1 - sin(a * xx))) * sin(a * xx) + kpare * (0.2e1 / 0.3e1 * (0.14e1 - sin(a * xx)) / Mref) ^ epn * cos(a * xx) * a / Mref / 0.45e2;
        end
end