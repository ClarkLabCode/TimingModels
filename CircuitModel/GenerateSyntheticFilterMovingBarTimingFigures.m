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

%% Generate temporal filters for plotting

% Second order exponential lowpass and highpass filters
[ params ] = SetModelParameters(baseArgIn{:});
t = (0:params.dt:1)';
tau = [0.075, 0.150, 0.225];
lp = 2 .* sqrt(params.dt) .*(tau.^(-3/2)) .* double(t>=0) .* t .* exp(-t./tau);
hp = 2 .* sqrt(params.dt) .* (tau.^(-3/2)) .* double(t>=0) .* (tau - t) .* exp(-t./tau);

MakeFigure;
hold on;
set(gca, 'colororder', corder([2,1,3],:));
plot(t, lp, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(num2str(1000*tau', '\\tau = %d ms'));
axis('square');
ConfAxis(16);
title('lowpass filters');

MakeFigure;
hold on;
set(gca, 'colororder', corder([2,1,3],:));
plot(t, hp, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(num2str(1000*tau', '\\tau = %d ms'));
axis('square');
ConfAxis(16);
title('highpass filters');

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% Center input time constant manipulations

% Set labels
legendStr = {'\tau_2 = 150 ms', '\tau_2 = 75 ms', '\tau_2 = 225 ms'};
titleStr = 'Model with single PD-offset delay';
labelStr = '';

% Standard model (tau2 = 150 ms)
[ ~, filters,  meanRespAllTau2_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.f2;

% Model with tau2 = 75 ms
[ ~, filters,  meanRespAllTau2_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau2', 0.075);
f(:,2) = filters.f2;

% Model with tau2 = 225 ms
[ ~, filters,  meanRespAllTau2_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau2', 0.225);
f(:,3) = filters.f2;

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
meanRespSel = cat(3, meanRespAllTau2_150(:,1:2), meanRespAllTau2_75(:,1:2), meanRespAllTau2_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% ND-offset delay input time constant manipulations with a LP filter

% Set labels
legendStr = {'\tau_1 = 150 ms', '\tau_1 = 75 ms', '\tau_1 = 225 ms'};
titleStr = 'Model with single ND-offset delay';
labelStr = '';

% Standard model (tau1 = 150 ms)
[ ~, filters,  meanRespAllTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.f1;

% Model with tau1 = 75 ms
[ ~, filters,  meanRespAllTau1_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.075);
f(:,2) = filters.f1;

% Model with tau1 = 150 ms
[ ~, filters,  meanRespAllTau1_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.225);
f(:,3) = filters.f1;

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
meanRespSel = cat(3, meanRespAllTau1_150(:,1:2), meanRespAllTau1_75(:,1:2), meanRespAllTau1_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% ND-offset delay input time constant manipulations with a LP filter and a parallel delay 

% Set labels
legendStr = {'\tau_1 = 150 ms', '\tau_1 = 75 ms', '\tau_1 = 225 ms'};
titleStr = 'Model with parallel ND-offset delays';
labelStr = '';

% Standard model (tau1 = 150 ms)
[ ~, filters,  meanRespAllTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'useParallelPdOffsetDelay', true);
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.f1;

% Model with tau1 = 75 ms
[ ~, filters,  meanRespAllTau1_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.075, 'useParallelPdOffsetDelay', true);
f(:,2) = filters.f1;

% Model with tau1 = 150 ms
[ ~, filters,  meanRespAllTau1_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.225, 'useParallelPdOffsetDelay', true);
f(:,3) = filters.f1;

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
meanRespSel = cat(3, meanRespAllTau1_150(:,1:2), meanRespAllTau1_75(:,1:2), meanRespAllTau1_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% ND-offset delay input time constant manipulations with a HP filter

% Set labels
legendStr = {'LP, \tau_1 = 150 ms','HP, \tau_1 = 150 ms', 'HP, \tau_1 = 75 ms', 'HP, \tau_1 = 225 ms'};
titleStr = 'Model with single ND-offset delay';
labelStr = '';

% Standard model (LP input 3, tau1 = 150 ms)
[ ~, filters,  meanRespAllLpTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),4);
f(:,1) = filters.f1;

% Model with HP input 1 and tau1 = 150 ms
[ ~, filters,  meanRespAllHpTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'hp1', true);
f(:,2) = filters.f1;

% Model with HP input 1 and tau1 = 75 ms
[ ~, filters,  meanRespAllHpTau1_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.075, 'hp1', true);
f(:,3) = filters.f1;

% Model with HP input 1 and tau1 = 150 ms
[ ~, filters,  meanRespAllHpTau1_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.225, 'hp1', true);
f(:,4) = filters.f1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllLpTau1_150(:,1:2), meanRespAllHpTau1_150(:,1:2), meanRespAllHpTau1_75(:,1:2), meanRespAllHpTau1_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4,2,3],:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% ND-offset delay input time constant manipulations with a HP filter and a parallel delay 

% Set labels
legendStr = {'LP, \tau_1 = 150 ms','HP, \tau_1 = 150 ms', 'HP, \tau_1 = 75 ms', 'HP, \tau_1 = 225 ms'};
titleStr = 'Model with parallel ND-offset delays';
labelStr = '';

% Standard model (LP input 3, tau1 = 150 ms)
[ ~, filters,  meanRespAllLpTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'useParallelPdOffsetDelay', true);
t = filters.t;
f = nan(length(t),4);
f(:,1) = filters.f1;

% Model with HP input 1 and tau1 = 150 ms
[ ~, filters,  meanRespAllHpTau1_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'hp1', true, 'useParallelPdOffsetDelay', true);
f(:,2) = filters.f1;

% Model with HP input 1 and tau1 = 75 ms
[ ~, filters,  meanRespAllHpTau1_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.075, 'hp1', true, 'useParallelPdOffsetDelay', true);
f(:,3) = filters.f1;

% Model with HP input 1 and tau1 = 150 ms
[ ~, filters,  meanRespAllHpTau1_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau1', 0.225, 'hp1', true, 'useParallelPdOffsetDelay', true);
f(:,4) = filters.f1;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllLpTau1_150(:,1:2), meanRespAllHpTau1_150(:,1:2), meanRespAllHpTau1_75(:,1:2), meanRespAllHpTau1_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4,2,3],:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;


%% PD-offset delay input time constant manipulations with a LP filter

% Set labels
legendStr = {'\tau_3 = 150 ms', '\tau_3 = 75 ms', '\tau_3 = 225 ms'};
titleStr = 'Model with single PD-offset delay';
labelStr = '';

% Standard model (tau3 = 150 ms)
[ ~, filters,  meanRespAllTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.f3;

% Model with tau3 = 75 ms
[ ~, filters,  meanRespAllTau3_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.075);
f(:,2) = filters.f3;

% Model with tau3 = 150 ms
[ ~, filters,  meanRespAllTau3_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.225);
f(:,3) = filters.f3;

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
meanRespSel = cat(3, meanRespAllTau3_150(:,1:2), meanRespAllTau3_75(:,1:2), meanRespAllTau3_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% PD-offset delay input time constant manipulations with a LP filter and a parallel delay 

% Set labels
legendStr = {'\tau_3 = 150 ms', '\tau_3 = 75 ms', '\tau_3 = 225 ms'};
titleStr = 'Model with parallel PD-offset delays';
labelStr = '';

% Standard model (tau3 = 150 ms)
[ ~, filters,  meanRespAllTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'useParallelPdOffsetDelay', true);
t = filters.t;
f = nan(length(t),3);
f(:,1) = filters.f3;

% Model with tau3 = 75 ms
[ ~, filters,  meanRespAllTau3_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.075, 'useParallelPdOffsetDelay', true);
f(:,2) = filters.f3;

% Model with tau3 = 150 ms
[ ~, filters,  meanRespAllTau3_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.225, 'useParallelPdOffsetDelay', true);
f(:,3) = filters.f3;

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
meanRespSel = cat(3, meanRespAllTau3_150(:,1:2), meanRespAllTau3_75(:,1:2), meanRespAllTau3_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder(1:3,:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% PD-offset delay input time constant manipulations with a HP filter

% Set labels
legendStr = {'LP, \tau_3 = 150 ms','HP, \tau_3 = 150 ms', 'HP, \tau_3 = 75 ms', 'HP, \tau_3 = 225 ms'};
titleStr = 'Model with single PD-offset delay';
labelStr = '';

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, filters,  meanRespAllLpTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:});
t = filters.t;
f = nan(length(t),4);
f(:,1) = filters.f3;

% Model with HP input 3 and tau3 = 150 ms
[ ~, filters,  meanRespAllHpTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'hp3', true);
f(:,2) = filters.f3;

% Model with HP input 3 and tau3 = 75 ms
[ ~, filters,  meanRespAllHpTau3_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.075, 'hp3', true);
f(:,3) = filters.f3;

% Model with HP input 3 and tau3 = 150 ms
[ ~, filters,  meanRespAllHpTau3_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.225, 'hp3', true);
f(:,4) = filters.f3;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllLpTau3_150(:,1:2), meanRespAllHpTau3_150(:,1:2), meanRespAllHpTau3_75(:,1:2), meanRespAllHpTau3_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4,2,3],:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

%% PD-offset delay input time constant manipulations with a HP filter and a parallel delay 

% Set labels
legendStr = {'LP, \tau_3 = 150 ms','HP, \tau_3 = 150 ms', 'HP, \tau_3 = 75 ms', 'HP, \tau_3 = 225 ms'};
titleStr = 'Model with parallel PD-offset delays';
labelStr = '';

% Standard model (LP input 3, tau3 = 150 ms)
[ ~, filters,  meanRespAllLpTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'useParallelPdOffsetDelay', true);
t = filters.t;
f = nan(length(t),4);
f(:,1) = filters.f3;

% Model with HP input 3 and tau3 = 150 ms
[ ~, filters,  meanRespAllHpTau3_150 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'hp3', true, 'useParallelPdOffsetDelay', true);
f(:,2) = filters.f3;

% Model with HP input 3 and tau3 = 75 ms
[ ~, filters,  meanRespAllHpTau3_75 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.075, 'hp3', true, 'useParallelPdOffsetDelay', true);
f(:,3) = filters.f3;

% Model with HP input 3 and tau3 = 150 ms
[ ~, filters,  meanRespAllHpTau3_225 ] = ThreeInputModelMovingBarResponses(v, config, barParam, baseArgIn{:}, 'tau3', 0.225, 'hp3', true, 'useParallelPdOffsetDelay', true);
f(:,4) = filters.f3;

% Plot filters for this model
MakeFigure;
hold on;
set(gca, 'colororder', corder([1,4,2,3],:));
plot(t, f, 'linewidth', 2);
plot(t, 0*t, '--k', 'linewidth', 2);
xlabel('time (s)');
ylabel('filter strength (contrast^{-1})');
legend(legendStr);
axis('square');
ConfAxis(16);
title(titleStr);

% Plot the results
meanRespSel = cat(3, meanRespAllLpTau3_150(:,1:2), meanRespAllHpTau3_150(:,1:2), meanRespAllHpTau3_75(:,1:2), meanRespAllHpTau3_225(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4,2,3],:)); 

% Clean up
clearvars -except config logV v barParam baseArgIn corder;

