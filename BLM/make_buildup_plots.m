%% figures to show build up as the period decreases to tau and below

clear Tlist tDelta fTrace;

tau = 100;
tauFast = 25;
t = [1:1e4];

[B,A]=butter(1,1/pi/tau,'low');
[BFast,AFast]= butter(1,1/pi/tauFast,'low');

Tlist = [400 200 100 50 25];

for ii=1:length(Tlist)
    tDelta = double(mod(t,Tlist(ii))==0);
    fTrace(ii,:) = filter(B,A,tDelta);
    fTraceFast(ii,:) = filter(BFast,AFast,tDelta);
end

figure; hold on;
cm=colormap('copper');
colororder(colorTriples(cm,length(Tlist)));
plot(t/tau-80,fTrace);
% colororder(colorTriples(cm,length(Tlist)));
% plot(t/tau-80,fTraceFast,':');
xlabel('time (units of tau)');
ylabel('response amplitude');
set(gca,'xlim',[0 12]);
title(['period/tau = ' num2str(Tlist/tau)]);
niceAxesLocal;

function out = colorTriples(initialMap,num)

n = size(initialMap,1);
out = interp1([0:(n-1)]/(n-1),initialMap,[0:(num-1)]/(num-1));

end

function niceAxesLocal

set(gca,'fontname','times','fontsize',20);
ax=findall(gca,'Type','line');
for i=1:length(ax)
    set(ax(i),'linewidth',2);
end
ax=findall(gcf,'Type','text');
for i=1:length(ax)
    set(ax(i),'fontname','times','fontsize',20);
end

end