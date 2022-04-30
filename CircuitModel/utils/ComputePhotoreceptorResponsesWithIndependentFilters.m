function [ resp1, resp2, resp3, blurArray ] = ComputePhotoreceptorResponsesWithIndependentFilters(stimArray, p, f)

% Compute the shift
prShift = floor(p.photoreceptorSpacing / p.dx);

% Blur the filter in space (assumes periodic boundary conditions)
if p.useSpatialFilter
    blurArray = fftshift(ifft(fft(f.spatialFilter,[],2) .* fft(stimArray,[],2), [], 2),2);
else
    blurArray = stimArray;
end

% Filter the spatially-blurred stimulus in time, using separate filters for
% each channel. Then, shift the stimulus to get each of the photoreceptor
% inputs. 
resp1 = circshift(filter(f.f1, 1, blurArray, [], 1), +prShift, 2);
resp2 = filter(f.f2, 1, blurArray, [], 1);
resp3 = circshift(filter(f.f3, 1, blurArray, [], 1), -prShift, 2);

end

