%% GENERATEMEASUREDFILTERMOVINGBARTIMINGFIGURES.M

addpath analyses models stimuli utils;

%% Set overall parameters

[ config ] = SetConfiguration('regenerateData',true);

baseArgIn = {
    'useMeasuredFilters', true,...   
    'deconvolveMeasuredFilters', false,...
    'smoothMeasuredFilters', true,...
    'filterSmoothingMethod', 'laguerre',...
    'filterLaguerreBasisOrder', 5,...
    'filterLaguerreBasisAlpha', 0.2,...
    'interpMethod', 'pchip'
    };

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

corder = [0,0,0; 0.915294117647059   0.281568627450980   0.287843137254902; 0.346666666666667   0.536000000000000   0.690666666666667;0.441568627450980   0.749019607843137   0.432156862745098];

%% Mi1 manipulations in model with Mi9, Mi1, Tm3, and Mi4

% Set labels
legendStr = {'WT', 'Mi1 > slo','Mi1 > sloRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and Mi4'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.Mi1;

% Mi1-slo
[ ~, filters, meanRespAllMi1slo ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo', 'gMi4', 0.3, 'gCT1', 0);
f(:,2) = filters.Mi1;

% Mi1-sloRNAi
[ ~, filters, meanRespAllMi1sloRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
f(:,3) = filters.Mi1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi1slo(:,1:2), meanRespAllMi1sloRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Mi1 manipulations in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'Mi1 > slo','Mi1 > sloRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.Mi1;

% Mi1-slo
[ ~, filters, meanRespAllMi1slo ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo');
f(:,2) = filters.Mi1;

% Mi1-sloRNAi
[ ~, filters, meanRespAllMi1sloRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi');
f(:,3) = filters.Mi1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi1slo(:,1:2), meanRespAllMi1sloRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Tm3 manipulations in model with Mi9, Mi1, Tm3, and Mi4

% Set labels
legendStr = {'WT', 'Tm3 > slo','Tm3 > sloRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and Mi4'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.Tm3;

% Tm3-slo
[ ~, filters, meanRespAllTm3slo ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo', 'gMi4', 0.3, 'gCT1', 0);
f(:,2) = filters.Tm3;

% Tm3-sloRNAi
[ ~, filters, meanRespAllTm3sloRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
f(:,3) = filters.Tm3;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllTm3slo(:,1:2), meanRespAllTm3sloRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Tm3 manipulations in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'Tm3 > slo','Tm3 > sloRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.Tm3;

% Tm3-slo
[ ~, filters, meanRespAllTm3slo ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo');
f(:,2) = filters.Tm3;

% Tm3-sloRNAi
[ ~, filters, meanRespAllTm3sloRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi');
f(:,3) = filters.Tm3;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllTm3slo(:,1:2), meanRespAllTm3sloRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Mi4 manipulations in model with Mi9, Mi1, Tm3, and Mi4

% Set labels
legendStr = {'WT', 'Mi4 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and Mi4'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.Mi4;

% Mi4-cacRNAi
[ ~, filters, meanRespAllMi4cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi', 'gMi4', 0.3, 'gCT1', 0);
f(:,2) = filters.Mi4;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi4cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Mi4 manipulations in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'Mi4 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.Mi4;

% Mi4-cacRNAi
[ ~, filters, meanRespAllMi4cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi');
f(:,2) = filters.Mi4;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi4cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% CT1 manipulations in model with Mi9, Mi1, Tm3, and CT1

% Set labels
legendStr = {'WT', 'CT1 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'gCT1', 0.3, 'gMi4', 0);
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.CT1;

% CT1-cacRNAi
[ ~, filters, meanRespAllCT1cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
f(:,2) = filters.CT1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllCT1cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% CT1 manipulations in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'CT1 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.CT1;

% CT1-cacRNAi
[ ~, filters, meanRespAllCT1cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi');
f(:,2) = filters.CT1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllCT1cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Mi9 manipulations in model with Mi9, Mi1, Tm3, and CT1

% Set labels
legendStr = {'WT', 'Mi9 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'gCT1', 0.3, 'gMi4', 0);
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.Mi9;

% Mi9-cacRNAi
[ ~, filters, meanRespAllCT1cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
f(:,2) = filters.Mi9;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllCT1cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Mi9 manipulations in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'Mi9 > cacRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, filters, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),2);
f(:,1) = filters.Mi9;

% Mi9-cacRNAi
[ ~, filters, meanRespAllCT1cacRNAi ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi');
f(:,2) = filters.Mi9;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllCT1cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;
