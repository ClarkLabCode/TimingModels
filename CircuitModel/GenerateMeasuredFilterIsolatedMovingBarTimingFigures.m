%% GENERATEMEASUREDFILTERISOLATEDMOVINGBARTIMINGFIGURES.M
% Stimulation with isolated bars, with compensation factors 

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
logV = (3:0.5:9)';
v = 2.^logV;

%% Configure plotting options

corder = [0,0,0; 0.915294117647059   0.281568627450980   0.287843137254902; 0.346666666666667   0.536000000000000   0.690666666666667;0.441568627450980   0.749019607843137   0.432156862745098];

%% Mi1 manipulations in model with Mi9, Mi1, Tm3, and Mi4

% Set labels
legendStr = {'WT', 'Mi1 > slo','Mi1 > sloRNAi'};
titleStr = {'Model with Mi9, Mi1, Tm3, and Mi4'};
labelStr = '';

meanRespAllWildType = nan(length(logV),4);
meanRespAllMi1slo = nan(length(logV),4);
meanRespAllMi1sloRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Mi1-slo
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo', 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllMi1slo(indV,:) = tempResp * 32 / v(indV);
        
        % Mi1-sloRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllMi1sloRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        
        % Mi1-slo
        [ ~, ~, meanRespAllMi1slo(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo', 'gMi4', 0.3, 'gCT1', 0);
        
        % Mi1-sloRNAi
        [ ~, ~, meanRespAllMi1sloRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllMi1slo = nan(length(logV),4);
meanRespAllMi1sloRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});

        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Mi1-slo
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo');
        meanRespAllMi1slo(indV,:) = tempResp * 32 / v(indV);
        
        % Mi1-sloRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi');
        meanRespAllMi1sloRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        
        % Mi1-slo
        [ ~, ~, meanRespAllMi1slo(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi1','Mi1_Slo');
        
        % Mi1-sloRNAi
        [ ~, ~, meanRespAllMi1sloRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameMi1','Mi1_SloRNAi');
    end
end


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

meanRespAllWildType = nan(length(logV),4);
meanRespAllTm3slo = nan(length(logV),4);
meanRespAllTm3sloRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Tm3-slo
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo', 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllTm3slo(indV,:) = tempResp * 32 / v(indV);
        
        % Tm3-sloRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllTm3sloRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        
        % Tm3-slo
        [ ~, ~, meanRespAllTm3slo(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo', 'gMi4', 0.3, 'gCT1', 0);
        
        % Tm3-sloRNAi
        [ ~, ~, meanRespAllTm3sloRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi', 'gMi4', 0.3, 'gCT1', 0);
    end
end

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
meanRespAllWildType = nan(length(logV),4);
meanRespAllTm3slo = nan(length(logV),4);
meanRespAllTm3sloRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Tm3-slo
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo');
        meanRespAllTm3slo(indV,:) = tempResp * 32 / v(indV);
        
        % Tm3-sloRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi');
        meanRespAllTm3sloRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        
        % Tm3-slo
        [ ~, ~, meanRespAllTm3slo(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameTm3','Tm3_Slo');
        
        % Tm3-sloRNAi
        [ ~, ~, meanRespAllTm3sloRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'fNameTm3','Tm3_SloRNAi');
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllMi4cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Mi4-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi', 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllMi4cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        
        % Mi4-cacRNAi
        [ ~, ~, meanRespAllMi4cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi', 'gMi4', 0.3, 'gCT1', 0);
        
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllMi4cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % Mi4-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi');
        meanRespAllMi4cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        
        % Mi4-cacRNAi
        [ ~, ~, meanRespAllMi4cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi4','Mi4_CacRNAi');
        
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllCT1cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % CT1-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
        meanRespAllCT1cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        
        % CT1-cacRNAi
        [ ~, ~, meanRespAllCT1cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
        
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllCT1cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % CT1-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi');
        meanRespAllCT1cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        
        % CT1-cacRNAi
        [ ~, ~, meanRespAllCT1cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameCT1','CT1_CacRNAi');
        
    end
end

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

meanRespAllWildType = nan(length(logV),4);
meanRespAllCT1cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % CT1-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
        meanRespAllCT1cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:}, 'gMi4', 0.3, 'gCT1', 0);
        
        % CT1-cacRNAi
        [ ~, ~, meanRespAllCT1cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi', 'gCT1', 0.3, 'gMi4', 0);
        
    end
end

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


meanRespAllWildType = nan(length(logV),4);
meanRespAllCT1cacRNAi = nan(length(logV),4);

for indV = 1:length(logV)
    if logV(indV) < 5
        
        barParam.barPeriod = 32;
        
        % WT
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        meanRespAllWildType(indV,:) = tempResp * 32 / v(indV);
        
        % CT1-cacRNAi
        [ ~, ~, tempResp ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi');
        meanRespAllCT1cacRNAi(indV,:) = tempResp * 32 / v(indV);
        
    else
        
        barParam.barPeriod = v(indV);
        
        % WT
        [ ~, ~, meanRespAllWildType(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:});
        
        % CT1-cacRNAi
        [ ~, ~, meanRespAllCT1cacRNAi(indV,:) ] = ThreeInputModelMovingBarResponses(v(indV), config, barParam, baseArgIn{:},'fNameMi9','Mi9_CacRNAi');
        
    end
end


% Plot the results
meanRespSel = cat(3, meanRespAllWildType(:,1:2), meanRespAllCT1cacRNAi(:,1:2));
PlotLogVelocityTuningCurvesWithLabels(logV,meanRespSel,legendStr,titleStr,labelStr,corder([1,4],:));

% Clean up
clearvars -except config logV v barParam baseArgIn corder;
