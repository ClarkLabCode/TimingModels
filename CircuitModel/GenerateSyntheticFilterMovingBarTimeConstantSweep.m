%% GENERATESYNTHETICFILTERMOVINGBARTIMINGFIGURES.M


addpath analyses models stimuli utils;

%% Set overall parameters

[ config ] = SetConfiguration('regenerateData', true);

baseArgIn = {'useMeasuredFilters', false};

addpath analyses models stimuli utils parameter-sweeps;

%% Set stimulus parameters

% Set bar parameters
barParam.barWidth = 5;
% barParam.barPeriod = 45;
barParam.barPeriod = 30;
barParam.mlum = 0;
barParam.c = 1;

% Define grid of velocities in log2-space
logV = (3:0.1:9)';
v = 2.^logV;

%% Configure plotting options

corder = [0,0,0;
    0.915294117647059,0.281568627450980,0.287843137254902;
    0.346666666666667,0.536000000000000,0.690666666666667;
    0.441568627450980,0.749019607843137,0.432156862745098
    ];

load('utils/blueRedColorMap.mat','cmpRed');

%% Set up parallel pool 


poolObj = gcp('nocreate');
if isempty(poolObj)
    poolObj = parpool('local');
end

%% Center input time constant manipulations

tau2 = (0.005:0.005:0.250)';

[ ~, ~,  meanRespWT ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});

meanRespAll = nan(length(v), length(tau2), 2);
parfor indT = 1:length(tau2)
    [ ~, ~,  tempResp ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau2', tau2(indT));
    meanRespAll(:,indT,:) = tempResp(:,1:2);
end

normPd = meanRespAll(:,:,1) ./ max(meanRespAll(:,:,1),[],1);
normNd = meanRespAll(:,:,2) ./ max(meanRespAll(:,:,1),[],1);

MakeFigure;
imagesc(tau2*1000, logV, normPd);
colormap(cmpRed);
hold on;
contour(tau2*1000, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('\tau_{2} (ms)');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'PD response normalized by max PD (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespWT(:,1) .* logV,1)) ./ squeeze(mean(meanRespWT(:,1),1));
plot(tau2*1000, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(tau2*1000, centerOfMassWt*ones(length(tau2),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);


MakeFigure;
imagesc(tau2*1000, logV, normNd);
colormap(cmpRed);
hold on;
contour(tau2*1000, logV, normNd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('\tau_{2} (ms)');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'ND response normalized by max PD (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normNd .* logV,1)) ./ squeeze(mean(normNd,1));
centerOfMassWt = squeeze(mean(meanRespWT(:,2) .* logV,1)) ./ squeeze(mean(meanRespWT(:,2),1));
plot(tau2*1000, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(tau2*1000, centerOfMassWt*ones(length(tau2),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);

%% Delayed input time constant manipulations with low-pass filter

tau3 = (0.005:0.005:0.250)';

[ ~, ~,  meanRespWT ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});

meanRespAll = nan(length(v), length(tau3), 2);
parfor indT = 1:length(tau3)
    [ ~, ~,  tempResp ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', tau3(indT));
    meanRespAll(:,indT,:) = tempResp(:,1:2);
end

normPd = meanRespAll(:,:,1) ./ max(meanRespAll(:,:,1),[],1);
normNd = meanRespAll(:,:,2) ./ max(meanRespAll(:,:,1),[],1);

MakeFigure;
imagesc(tau3*1000, logV, normPd);
colormap(cmpRed);
hold on;
contour(tau3*1000, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('\tau_{3} (ms)');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'PD response normalized by max PD (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespWT(:,1) .* logV,1)) ./ squeeze(mean(meanRespWT(:,1),1));
plot(tau3*1000, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(tau3*1000, centerOfMassWt*ones(length(tau3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model with single low-passed ND-offset delay');

MakeFigure;
imagesc(tau3*1000, logV, normNd);
colormap(cmpRed);
hold on;
contour(tau3*1000, logV, normNd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('\tau_{2} (ms)');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'ND response normalized by max PD (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normNd .* logV,1)) ./ squeeze(mean(normNd,1));
centerOfMassWt = squeeze(mean(meanRespWT(:,2) .* logV,1)) ./ squeeze(mean(meanRespWT(:,2),1));
plot(tau3*1000, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(tau3*1000, centerOfMassWt*ones(length(tau3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model with single low-passed ND-offset delay');