%% Set general parameters 
% Velocities
vRange = [-100:1:-1,1:1:100]'; 
% Spatial period, in units of detector separation
spatialPeriod = 6;
% Spatial filtering standard deviation
sigmaSpace = 1/sqrt(2*log(2))/2;
% Timestep (in ms)
dt = 0.01;
% Cutoff for periodic summation
M = 10;
% Amplitude of inhibitory arm
A = 2;
%% Plot model responses for varying excitatory time constants
% Time constants
tau1list = [40, 20, 60];
tau2list = [100, 100, 100];
% Plot model response
plotFirstOrderLowpassSmoothedBarlowLevickDiracCombResponse(vRange, spatialPeriod, dt, M, A, tau1list, tau2list, sigmaSpace);
%% Plot model responses for varying inhibitory time constants 
% Time constants
tau1list = [40, 40, 40];
tau2list = [100, 75, 125];
% Plot model response
plotFirstOrderLowpassSmoothedBarlowLevickDiracCombResponse(vRange, spatialPeriod, dt, M, A, tau1list, tau2list, sigmaSpace);
%% Utility functions 
function plotFirstOrderLowpassSmoothedBarlowLevickDiracCombResponse(vRange, spatialPeriod, dt, M, A, tau1list, tau2list, sigmaSpace)

if length(tau1list) ~= length(tau2list) 
    error('Time constant lists must be of equal length.');
end

% Iterate over time constants
tic;
meanResp = nan(length(vRange), length(tau1list));
for indT = 1:length(tau1list) 

    % Compute the response
    [meanResp(:,indT)] = computeBarlowLevickDiracCombResponseWithSpatialFilter(vRange, spatialPeriod, M, A, dt, tau1list(indT), tau2list(indT), sigmaSpace);

end
toc;

%% Plot the results

% corder = lines(length(tau1list));
corder = [0,0,0; 0.915294117647059   0.281568627450980   0.287843137254902; 0.346666666666667   0.536000000000000   0.690666666666667;0.441568627450980   0.749019607843137   0.432156862745098];
corder = corder(1:3,:);


MakeFigure;
colororder(corder);
plot(vRange, meanResp,'linewidth',2);
legend(strip(cellstr(num2str([tau1list',tau2list'],'\\tau_1 = %d ms, \\tau_2 = %d ms'))));
axis('square');
xlabel('velocity (units/s)');
ylabel('model response (arb. units)');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'Box','off');
legend('location', 'northwest');

MakeFigure;
colororder(corder);
plot(vRange, meanResp,'linewidth',2);
hold on;
plot(vRange, flipud(meanResp),'--','linewidth',2);
xlim([0, max(vRange)]);
legend(strip(cellstr(num2str([tau1list',tau2list'],'\\tau_1 = %d ms, \\tau_2 = %d ms'))));
axis('square');
xlabel('velocity (units/s)');
ylabel('model response (arb. units)');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'Box','off');
legend('location', 'northwest');

MakeFigure;
colororder(corder);
plot(vRange, meanResp - flipud(meanResp),'linewidth',2);
legend(strip(cellstr(num2str([tau1list',tau2list'],'\\tau_1 = %d ms, \\tau_2 = %d ms'))));
xlim([0, max(vRange)]);
axis('square');
xlabel('velocity (units/s)');
ylabel('PD - ND response (arb. units)');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'Box','off');

MakeFigure;
colororder(corder);
meanRespTest = meanResp ./ max(meanResp(vRange > 0,:),[],1);
meanRespTest(meanRespTest<0.001) = NaN;
% plot(vRange, meanResp ./ max(meanResp(vRange > 0,:),[],1),'linewidth',2);
plot(vRange, meanRespTest,'linewidth',2);
hold on;
% plot(vRange, flipud(meanResp) ./ max(meanResp(vRange > 0,:),[],1),'--','linewidth',2);
plot(vRange, flipud(meanRespTest),'--','linewidth',2);
xlim([0, max(vRange)]);
legend(strip(cellstr(num2str([tau1list',tau2list'],'\\tau_1 = %d ms, \\tau_2 = %d ms'))));
axis('square');
xlabel('velocity (units/s)');
ylabel('model response normalized by max PD');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'Box','off');
legend('location', 'northeast');


end

function [ z ] = gaussTail(x)
% Gaussian tail distribution function
z = erfc(x / sqrt(2)) / 2;
end

function [meanResp] = computeBarlowLevickDiracCombResponseWithSpatialFilter(vRange, spatialPeriod, M, A, dt, tau1, tau2, sigmaSpace)

% Allocate a container
meanResp = nan(length(vRange),1);

% Iterate over velocities
for indV = 1:length(vRange) 
    
    % Spacing between arms and smoothing scale, in ms
    deltaT = 1/vRange(indV)*1000;
    sV = sigmaSpace./abs(vRange(indV));

    % Period in time, assuming 30 degree spatial periods, in ms
    P = spatialPeriod/abs(vRange(indV))*1000;
    
    % One period in time, in ms
    tBase = (0:dt:P-dt)';

    % Define anonymous functions for filters 
    % (normalization is L1 in continuum measure, which is correct given how we average)
    f1 = @(t) gaussTail(sV/tau1 - t./sV).*exp(sV.^2/(2*tau1^2) - t/tau1)/tau1;
    f2 = @(t) gaussTail(sV/tau2 - t./sV).*exp(sV.^2/(2*tau2^2) - t/tau2)/tau2;
    
    % Generate the linear response via manual periodic summation
    % Note that this uses open boundary conditions
    linResp = zeros(length(tBase),1);
    for k = -M:M
        linResp = linResp + f1(tBase + k * P) - A * f2(tBase + k * P - deltaT);
    end

    % Compute the time-averaged nonlinear response
    meanResp(indV) = mean(max(0,linResp));
    
end

end

function MakeFigure()
figure('Position',[200,500,500,700],'WindowStyle','docked');
end