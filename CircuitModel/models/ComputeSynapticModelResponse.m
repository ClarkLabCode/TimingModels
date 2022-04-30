function [ meanResp, voltageResp, calciumResp ] = ComputeSynapticModelResponse(stimArray, params, filters)
%COMPUTESYNAPTICMODELRESPONSE

%% Define local utility functions

% Define the input rectifiers
if isinf(params.inputRectBeta)
    relu = @(z) (z .* (z>0));
else
    relu = @(z) z.*(erf(params.inputRectBeta * z)+1)/2;
end

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
resp1 = circshift(filter(filters.f1, 1, blurArray, [], 1), +prShift, 2);
resp2 = filter(filters.f2, 1, blurArray, [], 1);
resp3 = circshift(filter(filters.f3, 1, blurArray, [], 1), -prShift, 2);
if params.useParallelPdOffsetDelay
    resp4 = circshift(filter(filters.f4, 1, blurArray, [], 1), -prShift, 2);
end

%% Compute three-input conductance nonlinearity model response

% Compute each postsynaptic conductance
% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
% Input 4 ~ CT1
g1 = params.g1 .* relu(-resp1+params.threshold1);
g2 = params.g2 .* relu(+resp2-params.threshold2);
g3 = params.g3 .* relu(+resp3-params.threshold3);
if params.useParallelPdOffsetDelay
    g4 = params.g4 .* relu(+resp4-params.threshold4);
end

% Compute the numerator and denominator of the synaptic nonlinearity
if params.useParallelPdOffsetDelay
    numResp = params.V1 .* g1 + params.V2 .* g2 + params.V3 .* (g3 + g4)/2;
    denResp = params.gleak + g1 + g2 + (g3+g4)/2;
else
    numResp = params.V1 .* g1 + params.V2 .* g2 + params.V3 .* g3;
    denResp = params.gleak + g1 + g2 + g3;
end

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
calciumResp = halfsquare(voltageResp - params.outputThreshold);

%% Average model responses over time, space, and auxiliary phase 

meanResp = squeeze(nanmean(calciumResp(params.averagingMask,:,:,:),[4,2,1]));

end

