%%
t = 0:1e-3:10;
dt=0.001;
y1 = exp(-(t-0.045)/0.142);
y = (1-exp(-(t)/0.045))/(1-exp(-1));
kern = [y(t>=0 & t<=0.045) y1(t>=0.045 & t<2)];
figure; plot(dt:dt:numel(kern)*dt,kern)
x = zeros(size(t));
x(1:1000)=1;
y=conv(x,kern);
figure; plot(dt:dt:numel(y)*dt,y)