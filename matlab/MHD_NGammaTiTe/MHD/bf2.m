function f = bf2(X)
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
for i = 1:numel(x)
    xx = x(i);
    yy = y(i);
    f(:,1) = xx * (xx ^ 2 + yy ^ 2) / 0.15e2 + (0.1e1 + xx ^ 2 + yy ^ 2) * xx / 0.15e2 + yy * (xx ^ 2 + yy ^ 2) / 0.15e2 + (0.1e1 + xx ^ 2 + yy ^ 2) * yy / 0.15e2 - 0.899e3 / 0.225e3 * D;
    f(:,2) = xx * (xx ^ 2 + yy ^ 2) ^ 2 / 0.15e2 + 0.2e1 / 0.15e2 * (0.1e1 + xx ^ 2 + yy ^ 2) * (xx ^ 2 + yy ^ 2) * xx + yy * (xx ^ 2 + yy ^ 2) ^ 2 / 0.15e2 + 0.2e1 / 0.15e2 * (0.1e1 + xx ^ 2 + yy ^ 2) * (xx ^ 2 + yy ^ 2) * yy + Mref * (0.2e1 / 0.45e2 * xx * (0.10e2 + xx ^ 2 + yy ^ 2 - (xx ^ 2 + yy ^ 2) ^ 2 / 0.2e1) / Mref + (0.1e1 + xx ^ 2 + yy ^ 2) * (0.2e1 * xx - 0.2e1 * xx * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2 + 0.2e1 / 0.45e2 * xx * (0.10e2 + xx ^ 2 - yy ^ 2) / Mref + 0.2e1 / 0.45e2 * (0.1e1 + xx ^ 2 + yy ^ 2) * xx / Mref + 0.2e1 / 0.45e2 * yy * (0.10e2 + xx ^ 2 + yy ^ 2 - (xx ^ 2 + yy ^ 2) ^ 2 / 0.2e1) / Mref + (0.1e1 + xx ^ 2 + yy ^ 2) * (0.2e1 * yy - 0.2e1 * yy * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2 + 0.2e1 / 0.45e2 * yy * (0.10e2 + xx ^ 2 - yy ^ 2) / Mref - 0.2e1 / 0.45e2 * (0.1e1 + xx ^ 2 + yy ^ 2) * yy / Mref) - mu * (0.3596e4 / 0.225e3 * xx ^ 2 + 0.3596e4 / 0.225e3 * yy ^ 2 + 0.899e3 / 0.225e3 - 0.4e1 / 0.225e3 * yy * xx);
    f(:,3) = ((2 * xx * (10 + xx ^ 2 + yy ^ 2)) + (2 * (1 + xx ^ 2 + yy ^ 2) * xx) + 0.4e1 / 0.3e1 * xx * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) + 0.2e1 / 0.3e1 * (1 + xx ^ 2 + yy ^ 2) * (2 * xx - 2 * xx * (xx ^ 2 + yy ^ 2))) * (xx ^ 2 + yy ^ 2) / 0.30e2 + (((1 + xx ^ 2 + yy ^ 2) * (10 + xx ^ 2 + yy ^ 2)) + 0.2e1 / 0.3e1 * (1 + xx ^ 2 + yy ^ 2) * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1)) * xx / 0.15e2 + ((2 * yy * (10 + xx ^ 2 + yy ^ 2)) + (2 * (1 + xx ^ 2 + yy ^ 2) * yy) + 0.4e1 / 0.3e1 * yy * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) + 0.2e1 / 0.3e1 * (1 + xx ^ 2 + yy ^ 2) * (2 * yy - 2 * yy * (xx ^ 2 + yy ^ 2))) * (xx ^ 2 + yy ^ 2) / 0.30e2 + (((1 + xx ^ 2 + yy ^ 2) * (10 + xx ^ 2 + yy ^ 2)) + 0.2e1 / 0.3e1 * (1 + xx ^ 2 + yy ^ 2) * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1)) * yy / 0.15e2 - csii * (0.9889e4 / 0.225e3 + 0.3596e4 / 0.225e3 * (xx ^ 2) + 0.3596e4 / 0.225e3 * (yy ^ 2) - 0.4e1 / 0.225e3 * yy * xx) - kpari * ((0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) ^ epn * epn * (2 * xx - 2 * xx * (xx ^ 2 + yy ^ 2)) / (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) * ((2 * xx - 2 * xx * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2 + (2 * yy - 2 * yy * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2) / 0.30e2 + (0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) ^ epn * ((2 - 6 * xx ^ 2 - 2 * yy ^ 2) / Mref / 0.45e2 - 0.4e1 / 0.45e2 * yy * xx / Mref) / 0.30e2 + (0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) ^ epn * epn * (2 * yy - 2 * yy * (xx ^ 2 + yy ^ 2)) / (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) * ((2 * xx - 2 * xx * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2 + (2 * yy - 2 * yy * (xx ^ 2 + yy ^ 2)) / Mref / 0.45e2) / 0.30e2 + (0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) ^ epn * (-0.4e1 / 0.45e2 * yy * xx / Mref + (2 - 2 * xx ^ 2 - 6 * yy ^ 2) / Mref / 0.45e2) / 0.30e2) + pcourr * Mref * (xx ^ 2 + yy ^ 2) * (0.2e1 / 0.45e2 * xx * (10 + xx ^ 2 - yy ^ 2) / Mref + 0.2e1 / 0.45e2 * (1 + xx ^ 2 + yy ^ 2) * xx / Mref + 0.2e1 / 0.45e2 * yy * (10 + xx ^ 2 - yy ^ 2) / Mref - 0.2e1 / 0.45e2 * (1 + xx ^ 2 + yy ^ 2) * yy / Mref) + 0.3e1 / 0.4e1 * ((1 + xx ^ 2 + yy ^ 2) ^ 2) / tie * (0.2e1 / 0.3e1 * (10 + xx ^ 2 - yy ^ 2) / Mref - 0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) * sqrt(0.2e1) * sqrt(0.3e1) * ((10 + xx ^ 2 - yy ^ 2) / Mref) ^ (-0.3e1 / 0.2e1);
    f(:,4) = (xx * (10 + xx ^ 2 - yy ^ 2) * (xx ^ 2 + yy ^ 2)) / 0.9e1 + ((1 + xx ^ 2 + yy ^ 2) * (xx ^ 2 + yy ^ 2) * xx) / 0.9e1 + ((1 + xx ^ 2 + yy ^ 2) * (10 + xx ^ 2 - yy ^ 2) * xx) / 0.9e1 + (yy * (10 + xx ^ 2 - yy ^ 2) * (xx ^ 2 + yy ^ 2)) / 0.9e1 - ((1 + xx ^ 2 + yy ^ 2) * (xx ^ 2 + yy ^ 2) * yy) / 0.9e1 + ((1 + xx ^ 2 + yy ^ 2) * (10 + xx ^ 2 - yy ^ 2) * yy) / 0.9e1 - csie * (0.1798e4 / 0.45e2 + 0.899e3 / 0.75e2 * (xx ^ 2) - 0.899e3 / 0.75e2 * (yy ^ 2)) - kpare * ((0.2e1 / 0.3e1 * (10 + xx ^ 2 - yy ^ 2) / Mref) ^ epn * epn * xx / (10 + xx ^ 2 - yy ^ 2) * (0.2e1 / 0.45e2 * xx / Mref - 0.2e1 / 0.45e2 * yy / Mref) / 0.15e2 - (0.2e1 / 0.3e1 * (10 + xx ^ 2 - yy ^ 2) / Mref) ^ epn * epn * yy / (10 + xx ^ 2 - yy ^ 2) * (0.2e1 / 0.45e2 * xx / Mref - 0.2e1 / 0.45e2 * yy / Mref) / 0.15e2) - pcourr * Mref * (xx ^ 2 + yy ^ 2) * (0.2e1 / 0.45e2 * xx * (10 + xx ^ 2 - yy ^ 2) / Mref + 0.2e1 / 0.45e2 * (1 + xx ^ 2 + yy ^ 2) * xx / Mref + 0.2e1 / 0.45e2 * yy * (10 + xx ^ 2 - yy ^ 2) / Mref - 0.2e1 / 0.45e2 * (1 + xx ^ 2 + yy ^ 2) * yy / Mref) - 0.3e1 / 0.4e1 * ((1 + xx ^ 2 + yy ^ 2) ^ 2) / tie * (0.2e1 / 0.3e1 * (10 + xx ^ 2 - yy ^ 2) / Mref - 0.2e1 / 0.3e1 * (0.10e2 + (xx ^ 2) + (yy ^ 2) - ((xx ^ 2 + yy ^ 2) ^ 2) / 0.2e1) / Mref) * sqrt(0.2e1) * sqrt(0.3e1) * (((10 + xx ^ 2 - yy ^ 2) / Mref) ^ (-0.3e1 / 0.2e1));
    
end
