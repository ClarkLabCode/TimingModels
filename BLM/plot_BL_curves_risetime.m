function plot_BL_curves_risetime(tau1, tau2, tau1ilist, tau2ilist, velList)

% tau1 and 2 are the time constants for excitatory or inhibitory arms in ms
% velList is list of velocities to plot, in units per second
% width is width of acceptance function; 0 for delta, 1/2 for regular
% smoothing

% note 22/02/08 -- this script was having some trouble with the small time
% constants and long integration times, so I'm trying to fix that here.

%% script to play with BL signals and different velocities

A = 6; % additional weighting term for inhibitory component

T = [-2:1e-4:2]'; % time in s for simulation

% excitatory input always starts at 0

delta = 1; % detectors are 1 unit apart
sigmaSpace = 1/sqrt(2*log(2))/2; % spatial filtering sigma with FWHM at 1

% keyboard;

tau1i = tau1ilist(1);
tau2i = tau2ilist(1);

MakeFigure; 
subplot(2,1,1); hold on;
vRange = sort(velList(velList>0));
for ii=1:length(vRange)
    [excInput, inhInput] = compExcInhInputsRisetime(T, vRange(ii), tau1, tau2, tau1i, tau2i, A);    
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0 0.75 0],length(vRange),ii));
end
title(['delta function spatial ' num2str(sort(velList)) ' taus ' num2str([tau1,tau2])]);
niceAxesLocal;
ylabel('PD inputs');
set(gca,'xlim',400*[-.6 1],'ylim',[-40 14]);

subplot(2,1,2); hold on;
vRange = -sort(-velList(velList<0));
for ii=1:length(vRange)
    [excInput, inhInput] = compExcInhInputsRisetime(T, vRange(ii), tau1, tau2, tau1i, tau2i, A);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0.4940 0.1840 0.5560],length(vRange),ii));
end
xlabel('time (ms)');
set(gca,'xlim',400*[-.6 1],'ylim',[-40 14]);
ylabel('ND inputs');
niceAxesLocal;


MakeFigure; 
subplot(2,1,1); hold on;
vRange = sort(velList(velList>0));
for ii=1:length(vRange)
    [excInput, inhInput] = compExcInhInputsRisetimeWithSpatialFilter(T, vRange(ii), tau1, tau2, tau1i, tau2i, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0 0.75 0],length(vRange),ii));    
    legEnt{ii} = num2str(vRange(ii));
end
ylabel('PD responses');
set(gca,'xlim',400*[-.6 1],'ylim',[-40 14]);
title(['spatial filtering ' num2str(sort(velList)) ' taus ' num2str([tau1,tau2])]);
niceAxesLocal;

subplot(2,1,2); hold on;
vRange = -sort(-velList(velList<0));
for ii=1:length(vRange)
    [excInput, inhInput] = compExcInhInputsRisetimeWithSpatialFilter(T, vRange(ii), tau1, tau2, tau1i, tau2i, A, sigmaSpace);
    plot(T*1000,excInput + inhInput,'color',cmGraded([0 0 0],[0.4940 0.1840 0.5560],length(vRange),ii));
    legEnt{ii} = num2str(vRange(ii));
end
xlabel('time (ms)');
ylabel('ND responses');
set(gca,'xlim',400*[-.6 1],'ylim',[-40 14]);
niceAxesLocal;

%% plot tuning curves for with and without spatial filtering for a sweep of the tau1s in the tau list...

vRange = [-100:1:-1,1:1:100]; % range of velocities

rOut = nan(length(tau1ilist),length(vRange));
rOutSpatial = nan(length(tau1ilist),length(vRange));

for kk=1:length(tau1ilist)
    tau1 = tau1; % in ms
    tau2 = tau2; % in ms
    tau1i = tau1ilist(kk);
    tau2i = tau2ilist(1);

    [ excInput, inhInput ] = compExcInhInputsRisetime(T, vRange, tau1, tau2, tau1i, tau2i, A);
    rOut(kk,:) = mean(rect(excInput + inhInput),1);
    
    [ excInput, inhInput ] = compExcInhInputsRisetimeWithSpatialFilter(T, vRange, tau1, tau2, tau1i, tau2i, A, sigmaSpace);
    rOutSpatial(kk,:) = mean(rect(excInput + inhInput),1);
    
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
title(['delta spatial responses tau1i ' num2str(tau1ilist)]);
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

rOut = nan(length(tau2ilist),length(vRange));
rOutSpatial = nan(length(tau2ilist),length(vRange));

for kk=1:length(tau2ilist)
    tau1 = tau1; % in ms
    tau2 = tau2; % in ms
    tau1i = tau1ilist(1);
    tau2i = tau2ilist(kk);
    
    [ excInput, inhInput ] = compExcInhInputsRisetime(T, vRange, tau1, tau2, tau1i, tau2i, A);
    rOut(kk,:) = mean(rect(excInput + inhInput),1);
    
    [ excInput, inhInput ] = compExcInhInputsRisetimeWithSpatialFilter(T, vRange, tau1, tau2, tau1i, tau2i, A, sigmaSpace);
    rOutSpatial(kk,:) = mean(rect(excInput + inhInput),1);
    
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
title(['delta spatial responses tau2i ' num2str(tau2ilist)]);
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
end

function [excInput, inhInput] = compExcInhInputsRisetime(T, v, tau1, tau2, tau1i, tau2i, A)

deltaT = 1./v;

% Note that cutoffs in exponents are required to avoid numerical instability
excInput = exp(-rect(T)/tau1).*heaviside(T).*(1-exp(-rect(T)/tau1i))/(tau1-tau1*tau1i/(tau1+tau1i));
inhInput = -A*exp(-rect(T-deltaT)/tau2).*heaviside(T-deltaT).*(1-exp(-rect(T-deltaT)/tau2i))/(tau2-tau2*tau2i/(tau2+tau2i));

end

function [excInput, inhInput] = compExcInhInputsRisetimeWithSpatialFilter(T, v, tau1, tau2, tau1i, tau2i, A, sigma)

sV = sigma./abs(v);
deltaT = 1./v;

u1 = 1./(1./tau1 + 1./tau1i);
u2 = 1./(1./tau2 + 1./tau2i);

excInput = gaussTail(sV/tau1 - T./sV).*exp(sV.^2/(2*tau1^2) - T/tau1) ...
    - gaussTail(sV/u1 - T./sV).*exp(sV.^2/(2*u1^2) - T/u1);
excInput = excInput / (tau1-tau1*tau1i/(tau1+tau1i));
inhInput = - A*gaussTail(sV/tau2 - (T-deltaT)./sV).*exp(sV.^2/(2*tau2^2) - (T-deltaT)/tau2)...
    + A*gaussTail(sV/u2 - (T-deltaT)./sV).*exp(sV.^2/(2*u2^2) - (T-deltaT)/u2);
inhInput = inhInput / (tau2 - tau2*tau2i/(tau2+tau2i));

end

function [ z ] = gaussTail(x)
% Gaussian tail distribution function
z = erfc(x / sqrt(2)) / 2;
end

function [ z ] = hExp(x,a)
% Regulated version of (1/2).*erfc(x/sqrt(2)) .* exp(a.*x)



c = x < Inf; %
z = erfc( (x .* c) / sqrt(2)) .* exp(a .* (x .* c) ) .* c / 2;

end

function col3 = cmGraded(c1,c2,totNum,thisNum)
    
% returns graded color between c1 and c2 row triples
xi = [1,totNum];
col3 = interp1(xi,[c1(:)';c2(:)'],thisNum);
end

function x = foldDS(r)

x = r(:,end/2+1:end) - r(:,end/2:-1:1);
end

function x = addSpatialFilter(x,sigma,v) % convolves with temporal filter set to match spatial filter at velocity

% keyboard;
w = sigma/abs(v)*1000; % sigma in ms
t = [-round(w*3):round(w*3)]; % time in ms, out to 3 sigma
f = exp(-t.^2/w^2);
f = f/sum(f);
x = filtfilt(f,1,x); % two of these, each with w/sqrt(2), should get w back and has no phase shift
% keyboard;
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
