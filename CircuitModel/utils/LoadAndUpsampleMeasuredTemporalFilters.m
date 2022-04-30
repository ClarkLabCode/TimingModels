function [ f ] = LoadAndUpsampleMeasuredTemporalFilters(config, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Load filter data from file
load(fullfile(config.filterSourcePath), 'filterMat','dtFilter','filterList','tSec');

% Hard-coded filter timing
% dtFilter = 1/30;
f.raw.dt = dtFilter;

% Define meshes for interpolation
% WARNING: this will fail dramatically if the filters are not the same size
nT = length(tSec);
f.raw.t = tSec;
% f.t = (0:params.dt:dtFilter*(nT-1))';
f.t = (0:params.dt:1-params.dt)';

% Extract filters from structure
f.raw.fRawMi9 = filterMat(:,strcmp(filterList, params.fNameMi9));
f.raw.fRawMi1 = filterMat(:,strcmp(filterList, params.fNameMi1));
f.raw.fRawTm3 = filterMat(:,strcmp(filterList, params.fNameTm3));
f.raw.fRawMi4 = filterMat(:,strcmp(filterList, params.fNameMi4));
f.raw.fRawCT1 = filterMat(:,strcmp(filterList, params.fNameCT1));

% Deconvolve calcium kernels
if params.deconvolveMeasuredFilters
    
    % Concatenate filters together
    tempFilterMat = [f.raw.fRawMi9, f.raw.fRawMi1, f.raw.fRawTm3, f.raw.fRawMi4, f.raw.fRawCT1];
    
    % Deconvolve a first-order lowpass filter via division in FFT space
    lp = exp(-tSec./params.tauDec)./params.tauDec;
    ftdFilterMat = ifft(fft(tempFilterMat,[],1) ./ fft(lp,[],1),[],1);
    
    % Extract the results
    f.raw.fRawMi9 = ftdFilterMat(:,1);
    f.raw.fRawMi1 = ftdFilterMat(:,2);
    f.raw.fRawTm3 = ftdFilterMat(:,3);
    f.raw.fRawMi4 = ftdFilterMat(:,4);
    f.raw.fRawCT1 = ftdFilterMat(:,5);
    
end

% Smooth filters
if params.smoothMeasuredFilters
    
    switch params.filterSmoothingMethod
            
        case 'laguerre'
            
            % Get the desired subset of the polynomial basis
            laguerreFuncs = getLaguerrePolys(nT,params.filterLaguerreBasisOrder,...
                params.filterLaguerreBasisAlpha);
            laguerreProj = laguerreFuncs * laguerreFuncs';
            
            % Project the filters into the subspace
            f.raw.fSmoothMi9 = laguerreProj * f.raw.fRawMi9;
            f.raw.fSmoothMi1 = laguerreProj * f.raw.fRawMi1;
            f.raw.fSmoothTm3 = laguerreProj * f.raw.fRawTm3;
            f.raw.fSmoothMi4 = laguerreProj * f.raw.fRawMi4;
            f.raw.fSmoothCT1 = laguerreProj * f.raw.fRawCT1;
            
        otherwise
            error('Invalid filterSmoothingMethod: %s', params.filterSmoothingMethod);
    end
    
else
    
    f.raw.fSmoothMi9 = f.raw.fRawMi9;
    f.raw.fSmoothMi1 = f.raw.fRawMi1;
    f.raw.fSmoothTm3 = f.raw.fRawTm3;
    f.raw.fSmoothMi4 = f.raw.fRawMi4;
    f.raw.fSmoothCT1 = f.raw.fRawCT1;
    
end

% Interpolate filters
fMi9 = interp1(f.raw.t, f.raw.fSmoothMi9, f.t, params.interpMethod);
fMi1 = interp1(f.raw.t, f.raw.fSmoothMi1, f.t, params.interpMethod);
fTm3 = interp1(f.raw.t, f.raw.fSmoothTm3, f.t, params.interpMethod);
fMi4 = interp1(f.raw.t, f.raw.fSmoothMi4, f.t, params.interpMethod);
fCT1 = interp1(f.raw.t, f.raw.fSmoothCT1, f.t, params.interpMethod);

% L2-normalize filters
f.Mi9 = fMi9 ./ vecnorm(fMi9,2);
f.Mi1 = fMi1 ./ vecnorm(fMi1,2);
f.Tm3 = fTm3 ./ vecnorm(fTm3,2);
f.Mi4 = fMi4 ./ vecnorm(fMi4,2);
f.CT1 = fCT1 ./ vecnorm(fCT1,2);

end

