

addpath analyses models stimuli utils;

%% Set overall parameters

[ config ] = SetConfiguration('regenerateData',true);

baseArgIn = {
    'useMeasuredFilters', true,...   
    'deconvolveMeasuredFilters', true,...
    'tauDec', 0.250,...
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


%% Configure plotting options

corder = [0,0,0; 0.915294117647059   0.281568627450980   0.287843137254902; 0.346666666666667   0.536000000000000   0.690666666666667;0.441568627450980   0.749019607843137   0.432156862745098];

%% Mi4 and CT1 silencing in model with Mi9, Mi1, Tm3, Mi4, and CT1

% Set labels
legendStr = {'WT', 'Mi4 silenced','CT1 silenced','Mi4 and CT1 silenced'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1'};
labelStr = '';

% WT
[ ~, ~, logV, v, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:});

% Mi4 silenced
[ ~, ~, ~, ~, meanRespAllMi4Silenced ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:},'gMi4',0);

% CT1 silenced
[ ~, ~, ~, ~, meanRespAllCT1Silenced ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:},'gCT1',0);

% Both Mi4 and CT1 silenced
[ ~, ~, ~, ~, meanRespAllMi4AndCT1Silenced ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:},'gMi4',0,'gCT1',0);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi4Silenced(:,1:2), meanRespAllCT1Silenced(:,1:2), meanRespAllMi4AndCT1Silenced(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1:4],:));

% Clean up
clearvars -except config barParam baseArgIn corder;

%% Mi4 and CT1 silencing in model with Mi9, Mi1, Tm3, Mi4, and CT1, with compensation

% Set labels
legendStr = {'WT', 'Mi4 silenced','CT1 silenced'};
titleStr = {'Model with Mi9, Mi1, Tm3, Mi4, and CT1, with compensatory changes'};
labelStr = '';

% WT
[ ~, ~, logV, v, meanRespAllWildType ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:});

% Mi4 silenced
[ ~, ~, ~, ~, meanRespAllMi4Silenced ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:},'gMi4',0,'gCT1',0.3);

% CT1 silenced
[ ~, ~, ~, ~, meanRespAllCT1Silenced ] = ThreeInputModelMovingBarResponses(config, barParam, baseArgIn{:},'gCT1',0,'gMi4',0.3);

% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllMi4Silenced(:,1:2), meanRespAllCT1Silenced(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,2,3],:));

% Clean up
clearvars -except config barParam baseArgIn corder;
