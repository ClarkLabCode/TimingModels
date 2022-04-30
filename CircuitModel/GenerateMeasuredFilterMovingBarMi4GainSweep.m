
addpath analyses models stimuli utils;


% Base value of g3
g3_0 = 0.3;
g3_rel = (0:0.25:4)';
g3 = g3_0 * g3_rel;

%% Set plotting options

corder = cbrewer('div','Spectral', length(g3));
cmp = cbrewer('div','Spectral', 2^16);
load('utils/blueRedColorMap.mat','cmpRed');


%% Set overall parameters

[ config ] = SetConfiguration('regenerateData', true);

baseArgIn = {
    'useMeasuredFilters', true,...    
    'smoothMeasuredFilters', true,...
    'filterSmoothingMethod', 'laguerre',...
    'filterLaguerreBasisOrder', 5,...
    'filterLaguerreBasisAlpha', 0.2,...
    'interpMethod', 'pchip',...
    'tauMeasuredFilters', 0.05
    };

addpath analyses models stimuli utils parameter-sweeps;

%% Set stimulus parameters

% Set bar parameters
barParam.barWidth = 5;
% barParam.barPeriod = 45;
barParam.barPeriod = 30;
barParam.mlum = 0;
barParam.c = 1;


%% PD-offset delay input gain manipulations without a parallel delay, WT Mi4 filter

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, ~, logV, ~, meanRespAll ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'gMi4', g3_0, 'gCT1', 0);

% Iterate over other gains
respSweep = nan(length(logV), length(g3), 4);
for gInd = 1:length(g3)
    [ ~, ~, ~, ~, respSweep(:,gInd,:) ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'gCT1', 0, 'gMi4', g3(gInd));
end

normPd = respSweep(:,:,1) ./ max(respSweep(:,:,1),[],1);
normNd = respSweep(:,:,2) ./ max(respSweep(:,:,1),[],1);

MakeFigure;
hold on;
set(gca, 'colororder', corder);
colormap(cmp);
plot(logV, normNd, '--','linewidth',2);
plot(logV, normPd, '-','linewidth', 2);
plot(logV, meanRespAll(:,1) ./ max(meanRespAll(:,1)), '-k','linewidth',2);
plot(logV, meanRespAll(:,2) ./ max(meanRespAll(:,1)), '--k','linewidth',2);
cbar = colorbar;
caxis([min(g3_rel),max(g3_rel)]);
ylabel(cbar, 'g_3/g_{3,0}');
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
axis('square');
ConfAxis(16);
title('Model without a parallel delay, WT');

MakeFigure;
imagesc(g3_rel, logV, normPd);
colormap(cmpRed);
hold on;
contour(g3_rel, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('g_3/g_{3,0}');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'response (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespAll(:,1) .* logV,1)) ./ squeeze(mean(meanRespAll(:,1),1));
plot(g3_rel, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(g3_rel, centerOfMassWt*ones(length(g3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model without a parallel delay, WT');

%% PD-offset delay input gain manipulations with a parallel delay, WT Mi4 filter

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, filters, logV, v, meanRespAll ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:});

% Iterate over other gains
respSweep = nan(length(logV), length(g3), 4);
for gInd = 1:length(g3)
    [ ~, ~, ~, ~, respSweep(:,gInd,:) ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'gMi4', g3(gInd));
end

normPd = respSweep(:,:,1) ./ max(respSweep(:,:,1),[],1);
normNd = respSweep(:,:,2) ./ max(respSweep(:,:,1),[],1);

MakeFigure;
hold on;
set(gca, 'colororder', corder);
colormap(cmp);
plot(logV, normNd, '--','linewidth',2);
plot(logV, normPd, '-','linewidth', 2);
plot(logV, meanRespAll(:,1) ./ max(meanRespAll(:,1)), '-k','linewidth',2);
plot(logV, meanRespAll(:,2) ./ max(meanRespAll(:,1)), '--k','linewidth',2);
cbar = colorbar;
caxis([min(g3_rel),max(g3_rel)]);
ylabel(cbar, 'g_3/g_{3,0}');
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
axis('square');
ConfAxis(16);
title('Model with a parallel delay, WT');

MakeFigure;
imagesc(g3_rel, logV, normPd);
colormap(cmpRed);
hold on;
contour(g3_rel, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('g_3/g_{3,0}');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'response (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespAll(:,1) .* logV,1)) ./ squeeze(mean(meanRespAll(:,1),1));
plot(g3_rel, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(g3_rel, centerOfMassWt*ones(length(g3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model with a parallel delay, WT');

%% PD-offset delay input gain manipulations without a parallel delay, Mi4 > CacRNAi filter

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, ~, logV, ~, meanRespAll ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'gMi4', g3_0, 'gCT1', 0);

% Iterate over other gains
respSweep = nan(length(logV), length(g3), 4);
for gInd = 1:length(g3)
    [ ~, ~, ~, ~, respSweep(:,gInd,:) ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'fNameMi4','Mi4_CacRNAi', 'gMi4', g3(gInd), 'gCT1', 0);
end

normPd = respSweep(:,:,1) ./ max(respSweep(:,:,1),[],1);
normNd = respSweep(:,:,2) ./ max(respSweep(:,:,1),[],1);

MakeFigure;
hold on;
set(gca, 'colororder', corder);
colormap(cmp);
plot(logV, normNd, '--','linewidth',2);
plot(logV, normPd, '-','linewidth', 2);
plot(logV, meanRespAll(:,1) ./ max(meanRespAll(:,1)), '-k','linewidth',2);
plot(logV, meanRespAll(:,2) ./ max(meanRespAll(:,1)), '--k','linewidth',2);
cbar = colorbar;
caxis([min(g3_rel),max(g3_rel)]);
ylabel(cbar, 'g_3/g_{3,0}');
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
axis('square');
ConfAxis(16);
title('Model without a parallel delay, Mi4 > CacRNAi');

MakeFigure;
imagesc(g3_rel, logV, normPd);
colormap(cmpRed);
hold on;
contour(g3_rel, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('g_3/g_{3,0}');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'response (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespAll(:,1) .* logV,1)) ./ squeeze(mean(meanRespAll(:,1),1));
plot(g3_rel, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(g3_rel, centerOfMassWt*ones(length(g3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model without a parallel delay, Mi4 > CacRNAi');

%% PD-offset delay input gain manipulations with a parallel delay, Mi4 > CacRNAi filter

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, filters, logV, v, meanRespAll ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:});

% Iterate over other gains
respSweep = nan(length(logV), length(g3), 4);
for gInd = 1:length(g3)
    [ ~, ~, ~, ~, respSweep(:,gInd,:) ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:}, 'fNameMi4','Mi4_CacRNAi', 'gMi4', g3(gInd));
end

normPd = respSweep(:,:,1) ./ max(respSweep(:,:,1),[],1);
normNd = respSweep(:,:,2) ./ max(respSweep(:,:,1),[],1);

MakeFigure;
hold on;
set(gca, 'colororder', corder);
colormap(cmp);
plot(logV, normNd, '--','linewidth',2);
plot(logV, normPd, '-','linewidth', 2);
plot(logV, meanRespAll(:,1) ./ max(meanRespAll(:,1)), '-k','linewidth',2);
plot(logV, meanRespAll(:,2) ./ max(meanRespAll(:,1)), '--k','linewidth',2);
cbar = colorbar;
caxis([min(g3_rel),max(g3_rel)]);
ylabel(cbar, 'g_3/g_{3,0}');
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
axis('square');
ConfAxis(16);
title('Model with a parallel delay, Mi4 > CacRNAi');

MakeFigure;
imagesc(g3_rel, logV, normPd);
colormap(cmpRed);
hold on;
contour(g3_rel, logV, normPd, 0:0.1:1,'linewidth',2,'EdgeColor','k');
xlabel('g_3/g_{3,0}');
yticks(min(logV):1:max(logV));
ylim([min(logV),max(logV)]);
yticklabels(2.^(min(logV):1:max(logV))');
ylabel('velocity (\circ/s)');
cbar = colorbar;
ylabel(cbar,'response (arb. units)');
axis('xy','square');
ConfAxis(16);

centerOfMass = squeeze(mean(normPd .* logV,1)) ./ squeeze(mean(normPd,1));
centerOfMassWt = squeeze(mean(meanRespAll(:,1) .* logV,1)) ./ squeeze(mean(meanRespAll(:,1),1));
plot(g3_rel, centerOfMass,'-o','linewidth',2,'Color',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880],'MarkerEdgeColor','None');
plot(g3_rel, centerOfMassWt*ones(length(g3),1),'--','linewidth',2,'Color',[0.9290    0.6940    0.1250]);
title('Model with a parallel delay, Mi4 > CacRNAi');
