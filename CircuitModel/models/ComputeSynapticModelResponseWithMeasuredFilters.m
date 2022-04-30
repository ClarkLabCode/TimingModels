function [ meanResp, voltageResp, calciumResp, meanNumResp, meanDenResp, meanNumRespLN ] = ComputeSynapticModelResponseWithMeasuredFilters(stimArray, params, filters)
%COMPUTESYNAPTICMODELRESPONSE

%% Define local utility functions

% Define the input rectifiers
if isinf(params.inputRectBeta)
    relu = @(z) (z .* (z>0));
else
%     relu = @(z) z.*(erf(params.inputRectBeta * z)+1)/2;
    relu = @(z) log(1+exp(params.inputRectBeta * z))/params.inputRectBeta;
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
    if length(unique([filters.spatialFilterMi9,filters.spatialFilterMi1,filters.spatialFilterTm3,filters.spatialFilterMi4,filters.spatialFilterCT1])) == 1
        blurArray = fftshift(ifft(fft(filters.spatialFilterMi9,[],2) .* fft(stimArray,[],2), [], 2),2);
        blurArrayMi9 = blurArray;
        blurArrayMi1 = blurArray;
        blurArrayTm3 = blurArray;
        blurArrayMi4 = blurArray;
        blurArrayCT1 = blurArray;
    else
        blurArrayMi9 = fftshift(ifft(fft(filters.spatialFilterMi9,[],2) .* fft(stimArray,[],2), [], 2),2);
        blurArrayMi1 = fftshift(ifft(fft(filters.spatialFilterMi1,[],2) .* fft(stimArray,[],2), [], 2),2);
        blurArrayTm3 = fftshift(ifft(fft(filters.spatialFilterTm3,[],2) .* fft(stimArray,[],2), [], 2),2);
        blurArrayMi4 = fftshift(ifft(fft(filters.spatialFilterMi4,[],2) .* fft(stimArray,[],2), [], 2),2);
        blurArrayCT1 = fftshift(ifft(fft(filters.spatialFilterCT1,[],2) .* fft(stimArray,[],2), [], 2),2);
    end
else
    blurArrayMi9 = stimArray;
    blurArrayMi1 = stimArray;
    blurArrayTm3 = stimArray;
    blurArrayMi4 = stimArray;
    blurArrayCT1 = stimArray;
end

% Filter the spatially-blurred stimulus in time, using separate filters for
% each channel. Then, shift the stimulus to get each of the photoreceptor
% inputs.
respMi9 = circshift(filter(filters.Mi9, 1, blurArrayMi9, [], 1), +prShift, 2);
respMi1 = filter(filters.Mi1, 1, blurArrayMi1, [], 1);
respTm3 = filter(filters.Tm3, 1, blurArrayTm3, [], 1);
respMi4 = circshift(filter(filters.Mi4, 1, blurArrayMi4, [], 1), -prShift, 2);
respCT1 = circshift(filter(filters.CT1, 1, blurArrayCT1, [], 1), -prShift, 2);

%% Compute three-input conductance nonlinearity model response

% Compute each postsynaptic conductance
% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
gMi9 = params.gMi9 .* relu(-respMi9+params.threshold1);
gMi1 = params.gMi1 .* relu(+respMi1-params.threshold2);
gTm3 = params.gTm3 .* relu(+respTm3-params.threshold2);
gMi4 = params.gMi4 .* relu(+respMi4-params.threshold3);
gCT1 = params.gCT1 .* relu(+respCT1-params.threshold3);

% Compute the numerator and denominator of the three-input model
numResp = params.V1 .* gMi9 + params.V2 .* (gMi1 + gTm3) + params.V3 .* (gMi4+gCT1);
denResp = params.gleak + gMi9 + gMi1 + gTm3 + gMi4 + gCT1;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
calciumResp = halfsquare(voltageResp - params.outputThreshold);

%% Compute averaged numerator and denominator LNLN responses, if desired

% Factorization into a product of LNLN models
if nargout > 3
    meanNumResp = squeeze(nanmean(halfsquare(numResp(params.averagingMask,:,:,:)),[4,2,1]));
    meanDenResp = squeeze(nanmean(1./halfsquare(denResp(params.averagingMask,:,:,:)),[4,2,1]));
end

% Numerator LN model without intermediate rectification
if nargout > 5
    numRespLN = -params.V1*params.gMi9*respMi9 ...
        + params.V2 * (params.gMi1 * respMi1 + params.gTm3 * respTm3)...
        + params.V3 * (params.gMi4 * respMi4 + params.gCT1 * respCT1);
    meanNumRespLN = squeeze(nanmean(nanmean(nanmean(halfsquare(numRespLN(params.averagingMask,:,:,:)),4),2),1));
end

%% Average model responses over time and phase

meanResp = squeeze(nanmean(calciumResp(params.averagingMask,:,:,:),[4,2,1]));

end

