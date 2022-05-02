function plot_BL_curves(tau1list, tau2list, velList)

% tau1 and 2 are the time constants for excitatory or inhibitory arms in ms
% velList is list of velocities to plot, in units per second
% width is width of acceptance function; 0 for delta, 1/2 for regular
% smoothing

%% script to play with BL signals and different velocities

A = 6; % additional weighting term for inhibitory component

T = [-2:1e-4:2]'; % time in s for simulation

% excitatory input always starts at 0

delta = 1; % detectors are 1 unit apart
sigmaSpace = 1/sqrt(2*log(2))/2; % spatial filtering sigma with FWHM at 1

% keyboard;

tau1 = tau1list(1);
tau2 = tau2list(1); % lists are for the tuning curves below; plots are for first tau listed only...

MakeFigure; 
subplot(2,1,1); hold on;
vRange = sort(velList(velList>0));
for ii=1:length(vRange)
    [ excInput, inhInput ] = compExcInhInputs(T, vRange(ii), tau1, tau2, A);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0 0.75 0],length(vRange),ii));
end
title(['delta function spatial taus ' num2str([tau1,tau2])]);
niceAxesLocal;
ylabel('PD inputs');
legend(cellstr(num2str(round(vRange)', 'v = %-d')));
set(gca,'xlim',400*[-.6 1],'ylim',[-80 40]);

subplot(2,1,2); hold on;
vRange = -sort(-velList(velList<0));
for ii=1:length(vRange)
    [ excInput, inhInput ] = compExcInhInputs(T, vRange(ii), tau1, tau2, A);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0.4940 0.1840 0.5560],length(vRange),ii));
end
xlabel('time (ms)');
set(gca,'xlim',400*[-.6 1],'ylim',[-80 40]);
legend(cellstr(num2str(round(vRange)', 'v = %-d')));
ylabel('ND inputs');
niceAxesLocal;

MakeFigure; 
subplot(2,1,1); hold on;
vRange = sort(velList(velList>0));
for ii=1:length(vRange)
    [ excInput, inhInput ] = compExcInhInputsWithSpatialFilter(T, vRange(ii), tau1, tau2, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0 0.75 0],length(vRange),ii));
    legEnt{ii} = num2str(vRange(ii));
end
ylabel('PD responses');
set(gca,'xlim',400*[-.6 1],'ylim',[-40 14]);
title(['spatial filering taus ' num2str([tau1,tau2])]);
legend(cellstr(num2str(round(vRange)', 'v = %-d')));
niceAxesLocal;

subplot(2,1,2); hold on;
vRange = -sort(-velList(velList<0));
for ii=1:length(vRange)
    [ excInput, inhInput ] = compExcInhInputsWithSpatialFilter(T, vRange(ii), tau1, tau2, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color', cmGraded([0 0 0],[0.4940 0.1840 0.5560],length(vRange),ii));
    legEnt{ii} = num2str(vRange(ii));
end
xlabel('time (ms)');
ylabel('ND responses');
legend(cellstr(num2str(round(vRange)', 'v = %-d')));
set(gca,'xlim',400*[-.6 1],'ylim',[-60 14]);
niceAxesLocal;


%% plot tuning curves for with and without spatial filtering for a sweep of the tau1s in the tau list...
vRange = [-100:1:-1,1:1:100]; % range of velocities

rOut = nan(length(tau1list), length(vRange));
rOutSpatial = nan(length(tau1list), length(vRange));

for kk=1:length(tau1list)
    tau1 = tau1list(kk); % in ms
    tau2 = tau2list(1); % in ms

    [ excInput, inhInput ] = compExcInhInputs(T, vRange, tau1, tau2, A);
    rOut(kk,:) = mean(rect(excInput + inhInput), 1);
    
    [ excInput, inhInput ] = compExcInhInputsWithSpatialFilter(T, vRange, tau1, tau2, A, sigmaSpace);
    rOutSpatial(kk,:) = mean(rect(excInput + inhInput), 1);
    
    rOut(kk,:) = rOut(kk,:)/max(rOut(kk,:));
    rOutSpatial(kk,:) = rOutSpatial(kk,:)/max(rOutSpatial(kk,:));
end

% this falls off at high velocity becuase the inhibition starts to impinge
% on the excitation
MakeFigure; 
subplot(2,2,1); hold on;
colororder([0 0 0; 0.75 0 0; 0 0 0.75])
plot(vRange,rOut);
xlabel('velocity (units/s)');
ylabel('mean response');
title(['delta spatial responses tau1 ' num2str(tau1list)]);
niceAxesLocal;

subplot(2,2,2); hold on;
plot(vRange(end/2+1:end),foldDS(rOut)./max(foldDS(rOut),[],2));
xlabel('velocity (units/s)');
ylabel('mean response');
title('delta spatial PD - ND');
niceAxesLocal;

subplot(2,2,3); hold on;
plot(vRange,rOutSpatial);
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial responses');
niceAxesLocal;

subplot(2,2,4); hold on;
plot(vRange(end/2+1:end),foldDS(rOutSpatial)./max(foldDS(rOutSpatial),[],2));
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial PD - ND');
niceAxesLocal;

%% plot tuning curves for with and without spatial filtering for a sweep of the tau2s in the tau list...
clear rOut rOutSpatial;
vRange = [-100:1:-1,1:1:100]; % range of velocities

for kk=1:length(tau2list)
    tau1 = tau1list(1); % in ms
    tau2 = tau2list(kk); % in ms

    [ excInput, inhInput ] = compExcInhInputs(T, vRange, tau1, tau2, A);
    rOut(kk,:) = mean(rect(excInput + inhInput), 1);
    
    [ excInput, inhInput ] = compExcInhInputsWithSpatialFilter(T, vRange, tau1, tau2, A, sigmaSpace);
    rOutSpatial(kk,:) = mean(rect(excInput + inhInput), 1);
    
    rOut(kk,:) = rOut(kk,:)/max(rOut(kk,:));
    rOutSpatial(kk,:) = rOutSpatial(kk,:)/max(rOutSpatial(kk,:));
end

% this falls off at high velocity becuase the inhibition starts to impinge
% on the excitation
MakeFigure; 
colororder([0 0 0; 0.75 0 0; 0 0 0.75])
subplot(2,2,1); hold on;
plot(vRange,rOut);
xlabel('velocity (units/s)');
ylabel('mean response');
title(['delta spatial responses tau2 ' num2str(tau1list)]);
niceAxesLocal;

subplot(2,2,2); hold on;
plot(vRange(end/2+1:end),foldDS(rOut)./max(foldDS(rOut),[],2));
xlabel('velocity (units/s)');
ylabel('mean response');
title('delta spatial PD - ND');
niceAxesLocal;

subplot(2,2,3); hold on;
plot(vRange,rOutSpatial);
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial responses');
niceAxesLocal;

subplot(2,2,4); hold on;
plot(vRange(end/2+1:end),foldDS(rOutSpatial)./max(foldDS(rOutSpatial),[],2));
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial PD - ND');
niceAxesLocal;


%% now, just for kicks, let's do the same thing with delta functions (i.e., "the textbook model" we use as a straw man)

% excitatory input always starts at 0

delta = 1; % detectors are 1 unit apart
tau1 = 0;
tau2 = tau2list(1)-tau1list(1);
vRange = velList; % in units per second

MakeFigure;
vRange = sort(velList(velList>0));
subplot(2,1,1); hold on;
for ii=1:length(vRange)
    [ excInput, inhInput ] = compDiracTimeExcInhInputsWithSpatialFilter(T, vRange(ii), tau1, tau2, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0 0.75 0],length(vRange),ii));
end
xlabel('time (ms)');
ylabel('inputs');
niceAxesLocal;
title(['Gauss space, delta delay, vel = ' num2str(vRange)]);
set(gca,'xlim',400*[-.6 1]);



vRange = -sort(-velList(velList<0));
subplot(2,1,2); hold on;
for ii=1:length(vRange)
    [ excInput, inhInput ] = compDiracTimeExcInhInputsWithSpatialFilter(T, vRange(ii), tau1, tau2, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0.4940 0.1840 0.5560],length(vRange),ii));
end
xlabel('time (ms)');
ylabel('inputs');
title(['Gauss space, delta delay, vel = ' num2str(vRange)]);
niceAxesLocal;
set(gca,'xlim',400*[-.6 1]);


%% plot tuning curve -- this stays high everywhere, but falls low at single
% point when the inhibition and excitation cancel one another
vRange = [-100:1:-1,1:1:100]; % range of velocities
clear rOutSpatial;

rOutSpatial = nan(length(tau2list),length(vRange));

for kk=1:length(tau2list)
    tau1 = 0;
    tau2 = tau2list(kk)-tau1list(1); % just doing it this way
    
    [ excInput, inhInput ] = compDiracTimeExcInhInputsWithSpatialFilter(T, vRange, tau1, tau2, A, sigmaSpace);
    rOutSpatial(kk,:) = mean(rect(excInput + inhInput),1);
    
    rOutSpatial(kk,:) = rOutSpatial(kk,:)/max(rOutSpatial(kk,:));
end

MakeFigure; 
colororder([0 0 0; 0.75 0 0; 0 0 0.75])
subplot(2,1,1); hold on;
plot(vRange,rOutSpatial);
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial delta time, responses');
niceAxesLocal;legend('60 ms diff (100 ms inh)','35 ms diff (75 ms inh)','85 ms diff (125 ms inh)');

subplot(2,1,2); hold on;
plot(vRange(end/2+1:end),foldDS(rOutSpatial)./max(foldDS(rOutSpatial),[],2));
xlabel('velocity (units/s)');
ylabel('mean response');
title('+spatial delta time, PD - ND');
niceAxesLocal; legend('60 ms diff (100 ms inh)','35 ms diff (75 ms inh)','85 ms diff (125 ms inh)');

end

function [ excInput, inhInput ] = compExcInhInputs(T, v, tau1, tau2, A)

deltaT = 1./v; 
excInput = exp(-T/tau1).*heaviside(T)/tau1;
inhInput = -A*exp(-(T-deltaT)/tau2).*heaviside(T-deltaT)/tau2;
end

function [ excInput, inhInput ] = compExcInhInputsWithSpatialFilter(T, v, tau1, tau2, A, sigma)

if length(v)>1
%     keyboard;
end
sV = sigma./abs(v);
deltaT = 1./v;
excInput = gaussTail(sV/tau1 - T./sV).*exp(sV.^2/(2*tau1^2) - T/tau1)/tau1;
% excInput = excInput./sum(excInput);
if 1
    inhInput = - A*gaussTail(sV/tau2 - (T-deltaT)./sV).*exp(sV.^2/(2*tau2^2) - (T-deltaT)/tau2)/tau2;
else
    inhInput = gaussTail(sV/tau2 - (T-deltaT)./sV).*exp(sV.^2/(2*tau2^2) - (T-deltaT)/tau2)/tau2;
    inhInput = -A*inhInput./sum(inhInput);
end

end

function [ excInput, inhInput ] = compDiracTimeExcInhInputsWithSpatialFilter(T, v, tau1, tau2, A, sigma)

sV = sigma./abs(v);
deltaT = 1./v;

excInput = exp(-(T-tau1).^2 ./ (2.*sV.^2))./sqrt(2*pi*sV.^2);
inhInput = - A*exp(-(T-deltaT-tau2).^2 ./ (2.*sV.^2))./sqrt(2*pi*sV.^2);
end

function [ z ] = gaussTail(x)
% Gaussian tail distribution function
z = erfc(x / sqrt(2)) / 2;
end

function col3 = cmGraded(c1,c2,totNum,thisNum)
    
% returns graded color between c1 and c2 row triples
xi = [1,totNum];
col3 = interp1(xi,[c1(:)';c2(:)'],thisNum);
    
end

function x = foldDS(r)

x = r(:,end/2+1:end) - r(:,end/2:-1:1);

end

function x = rect(x)

x = max(0,x);

end

function MakeFigure()
figure('Position',[200,500,500,700],'WindowStyle','docked');
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