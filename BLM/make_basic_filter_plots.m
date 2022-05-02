%% makes some basic filter plots

t = [-50:500];

figure; hold on;
plot(t,exp(-t/20).*heaviside(t));
plot(t,exp(-t/40).*heaviside(t));
plot(t,exp(-t/60).*heaviside(t));
set(gca,'xlim',[-50 300]);
niceAxesLarge;

figure; hold on;
plot(t,exp(-t/75).*heaviside(t));
plot(t,exp(-t/100).*heaviside(t));
plot(t,exp(-t/125).*heaviside(t));
set(gca,'xlim',[-50 500]);
niceAxesLarge;

figure; hold on;
plot(t,exp(-t/100).*heaviside(t));
plot(t,exp(-t/40).*heaviside(t));
set(gca,'xlim',[-50 500]);
niceAxesLarge;

figure; hold on;
plot(t,t==75);
plot(t,t==100);
plot(t,t==125);
set(gca,'xlim',[-50 200]);
niceAxesLarge;

figure; hold on;
plot(t,t==0);
set(gca,'xlim',[-50 200]);
niceAxesLarge;

figure; hold on;
x=exp(-t/40).*(1-exp(-t/20)).*heaviside(t)/(40-800/60);
plot(t,x/max(x));
x=exp(-t/40).*(1-exp(-t/30)).*heaviside(t)/(40-1200/70);
plot(t,x/max(x));
x=exp(-t/40).*(1-exp(-t/40)).*heaviside(t)/(40-1600/80);
plot(t,x/max(x));
set(gca,'xlim',[-50 300]);
niceAxesLarge;

figure; hold on;
x = exp(-t/100).*(1-exp(-t/20)).*heaviside(t)/(100-2000/120);
plot(t,x/max(x));
x = exp(-t/100).*(1-exp(-t/40)).*heaviside(t)/(100-4000/140);
plot(t,x/max(x));
x = exp(-t/100).*(1-exp(-t/60)).*heaviside(t)/(100-6000/160);
plot(t,x/max(x));
set(gca,'xlim',[-50 500]);
niceAxesLarge;

figure; hold on;
x = exp(-t/100).*(1-exp(-t/40)).*heaviside(t)/(100-4000/140);
plot(t,x/max(x));
x=exp(-t/40).*(1-exp(-t/30)).*heaviside(t)/(40-1200/70);
plot(t,x/max(x));
set(gca,'xlim',[-50 500]);
niceAxesLarge;

