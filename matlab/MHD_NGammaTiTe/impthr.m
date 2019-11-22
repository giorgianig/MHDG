clear
close all
x = -10:0.001:10;
a = 1; 
delta = 0.1;
% f =a+(0.5+atan((x-a)/delta)/pi).*(x-a);
f = -delta*log(1./(1+exp( (x-a)/delta))) + a;
aux = (0.5+atan((x-a)/delta)/pi);

f1 = x;

df = exp((x-a)/delta)./(1+(exp( (x-a)/delta)));
figure, plot(x,f,'b',x,f1,'r--',x,a*ones(size(x)),'k--'); grid on
figure, plot(x,df); grid on