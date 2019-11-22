close all
% a = -1;

a =2;
x = 0:0.01:1;
y = (exp(a*x)-1)/(exp(a)-1);

% l = 1;
% x0 = 1;
% y = ( (1-(x-x0))/l).^a;
figure, plot(x,y),grid on
xlim([0,1])
ylim([0,1])
axis equal