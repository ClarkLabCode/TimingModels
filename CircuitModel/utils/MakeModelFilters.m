function [ f ] = MakeModelFilters(config, params)

%% Define temporal filters
% Note that these filters have unit L2 norm

if params.useMeasuredFilters
    
    % Load and upsample filters
    [ f ] = LoadAndUpsampleMeasuredTemporalFilters(config, params);
    
else
    
    % Truncate temporal vector to positive times
    t = params.t;
    t = t(t>=0);
    f.t = t;
    
    % Input 1 ~ Mi9
    if params.hp1
        f.f1 = sqrt(params.dt) .* 2 .* (params.tau1^(-3/2)) .* double(t>=0) .* (params.tau1 - t) .* exp(-t/params.tau1);
    else
        f.f1 = sqrt(params.dt) .* 2 .* (params.tau1^(-3/2)) .* double(t>=0) .* t .* exp(-t/params.tau1);
    end
    
    % Input 2 ~ Mi1/Tm3
    if params.hp2
        f.f2 = sqrt(params.dt) .* 2 .* (params.tau2^(-3/2)) .* double(t>=0) .* (params.tau2 - t) .* exp(-t/params.tau2);
    else
        f.f2 = sqrt(params.dt) .* 2 .* (params.tau2^(-3/2)) .* double(t>=0) .* t .* exp(-t/params.tau2);
    end
    
    % Input 3 ~ Mi4
    if params.hp3
        f.f3 = sqrt(params.dt) .* 2 .* (params.tau3^(-3/2)) .* double(t>=0) .* (params.tau3 - t) .* exp(-t/params.tau3);
    else
        f.f3 = sqrt(params.dt) .* 2 .* (params.tau3^(-3/2)) .* double(t>=0) .* t .* exp(-t/params.tau3);
    end
    
    % Input 4 ~ CT1
    if params.hp4
        f.f4 = sqrt(params.dt) .* 2 .* (params.tau4^(-3/2)) .* double(t>=0) .* (params.tau4 - t) .* exp(-t/params.tau4);
    else
        f.f4 = sqrt(params.dt) .* 2 .* (params.tau4^(-3/2)) .* double(t>=0) .* t .* exp(-t/params.tau4);
    end
    
end

%% Define spatial filters
% Note that these filter have unit L1 norm

if params.useMeasuredFilters
    
    sBlur = params.fwhmBlurMi9 / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilterMi9 = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
    sBlur = params.fwhmBlurMi1 / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilterMi1 = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
    sBlur = params.fwhmBlurTm3 / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilterTm3 = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
    sBlur = params.fwhmBlurMi4 / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilterMi4 = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
    sBlur = params.fwhmBlurCT1 / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilterCT1 = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
else
    
    % Compute spatial filter
    sBlur = params.fwhmBlur / (2*sqrt(2*log(2))); % Convert FWHM to STD
    f.spatialFilter = exp(-(params.x - params.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * params.dx;
    
end

end