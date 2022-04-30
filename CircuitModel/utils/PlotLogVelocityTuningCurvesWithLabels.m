function PlotLogVelocityTuningCurvesWithLabels(logV,meanResp,legendStr,titleStr,labelStr,corder)

% Validate inputs
if size(meanResp,3) ~= size(corder,1)
    error('One color must be specified for each tuning curve!');
end

%% Plot PD and ND responses

% Plot raw responses
MakeFigure;
hold on;
for ind = 1:size(meanResp,3) % PD
    plot(logV, meanResp(:,1,ind), '-','lineWidth', 2, 'color', corder(ind,:));
end
for ind = 1:size(meanResp,3) % ND
    plot(logV, meanResp(:,2,ind), '--','lineWidth', 2, 'color', corder(ind,:));
end

legend(legendStr);
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
axis('square');
ConfAxis(16);

title(strcat(titleStr,labelStr));

% Plot responses normalized by max PD (of full model)
normResp = meanResp ./ max(meanResp(:,1,:),[],1);
MakeFigure;
hold on;
for ind = 1:size(meanResp,3) % PD
    plot(logV, normResp(:,1,ind), '-','lineWidth', 2, 'color', corder(ind,:));
end
for ind = 1:size(meanResp,3) % ND
    plot(logV, normResp(:,2,ind), '--','lineWidth', 2, 'color', corder(ind,:));
end

legend(legendStr);
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('response normalized by max PD');
axis('square');
ConfAxis(16);
title(strcat(titleStr,labelStr));

% Plot centers of mass in log-velocity space
centerOfMass = squeeze(mean(meanResp(:,1,:) .* logV,1)) ./ squeeze(mean(meanResp(:,1,:),1));

MakeFigure;
hold on;
bar(1:size(meanResp,3), centerOfMass, 'FaceColor', 'flat','EdgeColor','none', 'CData', corder);
xticks(1:size(meanResp,3));
xticklabels(legendStr);
yticks(log2([16;20;24;28;32]));
ylim([4 5]);
yticklabels([16;20;24;28;32]);
ylabel('log-velocity center of mass');
axis('square');
ConfAxis(16);

title(strcat(titleStr,labelStr));

%% Plot PD - ND responses

% Plot raw responses
MakeFigure;
hold on;
for ind = 1:size(meanResp,3) 
    plot(logV, meanResp(:,1,ind) - meanResp(:,2,ind), '-','lineWidth', 2, 'color', corder(ind,:));
end

legend(legendStr);
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('PD - ND response (arb. units)');
axis('square');
ConfAxis(16);

title(strcat(titleStr,labelStr));

% Plot responses normalized by max PD - ND
normResp = (meanResp(:,1,:) - meanResp(:,2,:)) ./ max(meanResp(:,1,:) - meanResp(:,2,:),[],1);
MakeFigure;
hold on;
for ind = 1:size(meanResp,3) % PD
    plot(logV, normResp(:,1,ind), '-','lineWidth', 2, 'color', corder(ind,:));
end

legend(legendStr);
xticks(min(logV):1:max(logV));
xlim([min(logV),max(logV)]);
xticklabels(2.^(min(logV):1:max(logV))');
xlabel('velocity (\circ/s)');
ylabel('PD - ND response normalized by max PD - ND');
axis('square');
ConfAxis(16);
title(strcat(titleStr,labelStr));

end