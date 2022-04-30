function [ meanResp, voltageResp, calciumResp ] = ComputeLNModelResponseWithMeasuredFilters(stimArray, params, filters)
%COMPUTESYNAPTICMODELRESPONSE

%% Define local utility functions

% Define the output half-quadratic
if isinf(params.outputRectBeta)
    halfsquare = @(z) (z .* (z>0)).^2;
else
    halfsquare = @(z) (z.*(erf(params.outputRectBeta*z)+1)/2).^2;
end

%% Compute photoreceptor responses

% Compute the shift
prShift = floor(params.photoreceptorSpacing / params.dx);

% Blur the filter in space (assumes periodic boundary conditions)
if params.useSpatialFilter
    blurArray = fftshift(ifft(fft(filters.spatialFilter,[],2) .* fft(stimArray,[],2), [], 2),2);
else
    blurArray = stimArray;
end

% Filter the spatially-blurred stimulus in time, using separate filters for
% each channel. Then, shift the stimulus to get each of the photoreceptor
% inputs. 
respMi9 = circshift(filter(filters.Mi9, 1, blurArray, [], 1), +prShift, 2);
respMi1 = filter(filters.Mi1, 1, blurArray, [], 1);
respTm3 = filter(filters.Tm3, 1, blurArray, [], 1);
respMi4 = circshift(filter(filters.Mi4, 1, blurArray, [], 1), -prShift, 2);
respCT1 = circshift(filter(filters.CT1, 1, blurArray, [], 1), -prShift, 2);

%% Compute three-input conductance nonlinearity model response

% Compute each postsynaptic conductance
% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
gMi9 = params.gMi9 .* (-respMi9+params.threshold1);
gMi1 = params.gMi1 .* (+respMi1-params.threshold2);
gTm3 = params.gTm3 .* (+respTm3-params.threshold2);
gMi4 = params.gMi4 .* (+respMi4-params.threshold3);
gCT1 = params.gCT1 .* (+respCT1-params.threshold3);

% Compute the numerator and denominator of the three-input model
voltageResp = params.V1 .* gMi9 + params.V2 .* (gMi1 + gTm3) + params.V3 .* (gMi4+gCT1);

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
calciumResp = halfsquare(voltageResp - params.outputThreshold);

%% Average model responses over time and phase

meanResp = squeeze(nanmean(calciumResp(params.averagingMask,:,:,:),[4,2,1]));

end

