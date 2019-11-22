function f = bf5(X)
global testcase Mesh diff_n diff_u diff_ei diff_ee diff_pari diff_pare tie axisym neq  epn pcourr
global Mref
D = diff_n;
mu = diff_u;
csii = diff_ei;
csie = diff_ee;
kpari = diff_pari;
kpare = diff_pare;



% allocation
f = zeros(size(X,1),neq);
% assign
x = X(:,1);
y = X(:,2);




if axisym,error('This is NOT an axisymmetric case!'),end
a = pi/2;
for i = 1:numel(x)
    xx = x(i);
    yy = y(i);
    f(:,1) = -cos(a * xx) * a * sin(a * xx) / 0.30e2 + (0.11e1 - sin(a * xx)) * cos(a * xx) * a / 0.30e2 - 0.899e3 / 0.900e3 * D * sin(a * xx) * a ^ 2;
    f(:,2) = -cos(a * xx) * a * sin(a * xx) ^ 2 / 0.30e2 + (0.11e1 - sin(a * xx)) * sin(a * xx) * cos(a * xx) * a / 0.15e2 + Mref * (-cos(a * xx) * a * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref / 0.45e2 + (0.11e1 - sin(a * xx)) * (-cos(a * xx) * a - cos(a * xx) * a * sin(a * xx)) / Mref / 0.45e2 - cos(a * xx) * a * (0.14e1 - sin(a * xx)) / Mref / 0.45e2 - (0.11e1 - sin(a * xx)) * cos(a * xx) * a / Mref / 0.45e2) - mu * (0.899e3 / 0.900e3 * sin(a * xx) ^ 2 * a ^ 2 - 0.899e3 / 0.450e3 * cos(a * xx) ^ 2 * a ^ 2 - 0.899e3 / 0.900e3 * (0.11e1 - sin(a * xx)) * sin(a * xx) * a ^ 2);
    f(:,3) = (-cos(a * xx) * a * (0.2e1 - sin(a * xx)) - (0.11e1 - sin(a * xx)) * cos(a * xx) * a - 0.2e1 / 0.3e1 * cos(a * xx) * a * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) + 0.2e1 / 0.3e1 * (0.11e1 - sin(a * xx)) * (-cos(a * xx) * a - cos(a * xx) * a * sin(a * xx))) * sin(a * xx) / 0.30e2 + ((0.11e1 - sin(a * xx)) * (0.2e1 - sin(a * xx)) + 0.2e1 / 0.3e1 * (0.11e1 - sin(a * xx)) * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1)) * cos(a * xx) * a / 0.30e2 - csii * (0.899e3 / 0.900e3 * sin(a * xx) * a ^ 2 * (0.2e1 - sin(a * xx)) + 0.899e3 / 0.450e3 * cos(a * xx) ^ 2 * a ^ 2 + 0.899e3 / 0.900e3 * (0.11e1 - sin(a * xx)) * sin(a * xx) * a ^ 2) - kpari * ((0.2e1 / 0.3e1 * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref) ^ epn * epn * (-cos(a * xx) * a - cos(a * xx) * a * sin(a * xx)) ^ 2 / (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref / 0.1350e4 + (0.2e1 / 0.3e1 * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref) ^ epn * (sin(a * xx) * a ^ 2 + sin(a * xx) ^ 2 * a ^ 2 - cos(a * xx) ^ 2 * a ^ 2) / Mref / 0.1350e4) + pcourr * Mref * sin(a * xx) * (-cos(a * xx) * a * (0.14e1 - sin(a * xx)) / Mref / 0.45e2 - (0.11e1 - sin(a * xx)) * cos(a * xx) * a / Mref / 0.45e2) + 0.3e1 / 0.4e1 * (0.11e1 - sin(a * xx)) ^ 2 * (0.2e1 / 0.3e1 * (0.14e1 - sin(a * xx)) / Mref - 0.2e1 / 0.3e1 * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref) / tie * sqrt(0.2e1) * sqrt(0.3e1) * ((0.14e1 - sin(a * xx)) / Mref) ^ (-0.3e1 / 0.2e1);
    
    f(:,4) = -cos(a * xx) * a * (0.14e1 - sin(a * xx)) * sin(a * xx) / 0.18e2 - (0.11e1 - sin(a * xx)) * sin(a * xx) * cos(a * xx) * a / 0.18e2 + (0.11e1 - sin(a * xx)) * (0.14e1 - sin(a * xx)) * cos(a * xx) * a / 0.18e2 - csie * (0.899e3 / 0.900e3 * sin(a * xx) * a ^ 2 * (0.14e1 - sin(a * xx)) + 0.899e3 / 0.450e3 * cos(a * xx) ^ 2 * a ^ 2 + 0.899e3 / 0.900e3 * (0.11e1 - sin(a * xx)) * sin(a * xx) * a ^ 2) - kpare * ((0.2e1 / 0.3e1 * (0.14e1 - sin(a * xx)) / Mref) ^ epn * epn * cos(a * xx) ^ 2 * a ^ 2 / (0.14e1 - sin(a * xx)) / Mref / 0.1350e4 + (0.2e1 / 0.3e1 * (0.14e1 - sin(a * xx)) / Mref) ^ epn * sin(a * xx) * a ^ 2 / Mref / 0.1350e4) - pcourr * Mref * sin(a * xx) * (-cos(a * xx) * a * (0.14e1 - sin(a * xx)) / Mref / 0.45e2 - (0.11e1 - sin(a * xx)) * cos(a * xx) * a / Mref / 0.45e2) - 0.3e1 / 0.4e1 * (0.11e1 - sin(a * xx)) ^ 2 * (0.2e1 / 0.3e1 * (0.14e1 - sin(a * xx)) / Mref - 0.2e1 / 0.3e1 * (0.2e1 - sin(a * xx) - sin(a * xx) ^ 2 / 0.2e1) / Mref) / tie * sqrt(0.2e1) * sqrt(0.3e1) * ((0.14e1 - sin(a * xx)) / Mref) ^ (-0.3e1 / 0.2e1);
    
end
