
%%

corder = lines(2);

%% Set overall parameters

[ config ] = SetConfiguration('regenerateData',true);

baseArgIn = {
    'useMeasuredFilters', true,...   
    'deconvolveMeasuredFilters', false,...
    'smoothMeasuredFilters', true,...
    'filterSmoothingMethod', 'laguerre',...
    'filterLaguerreBasisOrder', 5,...
    'filterLaguerreBasisAlpha', 0.2,...
    'interpMethod', 'pchip',...
    'inputRectBeta', 1
    };

%%

load(config.filterSourcePath, 'filterMat','filterSem', 'filterLabel','filterList','tSec','dtFilter');

% Get the desired subset of the polynomial basis
laguerreFuncs = getLaguerrePolys(length(tSec), 5, 0.2);
laguerreProj = laguerreFuncs * laguerreFuncs';

lagFilterMat = laguerreProj * filterMat;

%%

for ind = 1:length(filterLabel)
    
    MakeFigure;
    hold on;
    plot(tSec, filterMat(:,ind) ./ sum(abs(filterMat(:,ind))) ./ dtFilter, 'linewidth', 2, 'Color', corder(1,:));
    plot(tSec, lagFilterMat(:,ind) ./ sum(abs(lagFilterMat(:,ind))) ./ dtFilter, 'linewidth', 2, 'Color', corder(2,:));
    plot(tSec, (filterMat(:,ind) + filterSem(:,ind)) ./ sum(abs(filterMat(:,ind))) ./ dtFilter, 'linewidth', 1, 'Color', corder(1,:));
    plot(tSec, (filterMat(:,ind) - filterSem(:,ind)) ./ sum(abs(filterMat(:,ind))) ./ dtFilter, 'linewidth', 1, 'Color', corder(1,:));
    
    axis('square');
    xlim([0 1]);
%     ylim([-1 1.5]);
    title(filterLabel{ind});
    ConfAxis(16);
    
end