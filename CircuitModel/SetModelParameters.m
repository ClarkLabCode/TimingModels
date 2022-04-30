function [ p ] = SetModelParameters(varargin)

%% Set whether to use measured temporal filters

% Global flag
p.useMeasuredFilters = false;

%% Set timing parameters

% Interleave duration (used to pad beginning and end of stimulus)
p.tInt = 0;

% Stimulus duration
% p.tOn = 3;
p.tOn = 5;

% Duration to exclude onset transients
% p.tAv = 1;
p.tAv = 0;

% Temporal resolution (s)
p.dt = 1/240;

%% Set spatial parameters

% Spatial extent (fills entire 360-degree space)
p.xExtent = 360;

p.xPad = 30;

% Spatial resolution (degrees)
p.dx = 0.5;

%% Set parameters for spatial filters

p.useSpatialFilter = true;

% Blur filter FWHM (degrees)
p.fwhmBlur = 5.7;

p.fwhmBlurMi9 = 5.7;
p.fwhmBlurMi1 = 5.7;
p.fwhmBlurTm3 = 5.7;
p.fwhmBlurMi4 = 5.7;
p.fwhmBlurCT1 = 5.7;

% Averaging filter FWHM (degrees)
p.fwhmAverage = 20 * (2*sqrt(2*log(2)));

% Set photoreceptor spacing (degrees)
p.photoreceptorSpacing = 5;
% p.photoreceptorSpacing = 5.1;

%% Set temporal filter parameters

% Filter time constants (s)
p.tau1 = 0.15;
p.tau2 = 0.15;
p.tau3 = 0.15;
p.tau4 = 0.15;

% Select whether to use 2nd order low- or high-pass filters
p.hp1 = false;
p.hp2 = true;
p.hp3 = false;
p.hp4 = false;

% Select whether to use a parallel PD-offset delay arm
p.useParallelPdOffsetDelay = false;

%% Set spatial phase shifts
% Note that these phase shifts are not used in the linearity analysis,
% where the phase shifts are fixed to be (0:1:7)/8*pi as in Wienecke et al.

p.useRandomShifts = false;

%% Set rectification parameters

% For soft rectification:
% params.inputRectBeta = 2;
% params.outputRectBeta = 32;

% For hard rectification:
p.inputRectBeta = Inf;
p.outputRectBeta = Inf;

% Input rectification thresholds
% Note that there is a sign-inversion for input 1
p.threshold1 = 0;
p.threshold2 = 0;
p.threshold3 = 0;
p.threshold4 = 0; % Only used if useParallelPdOffsetDelay is true

% Output rectification threshold
p.outputThreshold = 0;

%% Set reversal potentials & conductances

% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4

% Reversal potentials
p.V1 = - 30;
p.V2 = + 60;
p.V3 = - 30;

% Leak conductance
p.gleak = 1;

% Conductance gains 
p.g1 = 0.3;
p.g2 = 0.1;
p.g3 = 0.3;
p.g4 = 0.3; % Only used if useParallelPdOffsetDelay is true

%% Set parameters relating to usage of measured temporal filters

% Flag for whether or not to smooth filters
p.smoothMeasuredFilters = true;

% Flag for whether or not to deconvolve calcium kernels
p.deconvolveMeasuredFilters = true;

% Parameters for deconvolution 
% (note: defaults to a scalar, but can input a row vector in the order Mi9,
% Mi1, Tm3, Mi4, CT1 if different time constants are desired)
p.tauDec = 0.250;

% Parameters for Laguerre basis smoothing
p.filterSmoothingMethod = 'laguerre';
p.filterLaguerreBasisOrder = 5;
p.filterLaguerreBasisAlpha = 0.2;

% Method to interpolate to simulation temporal resolution
p.interpMethod = 'linear';

% Selectors
p.fNameMi9 = 'Mi9';
p.fNameMi1 = 'Mi1';
p.fNameMi4 = 'Mi4';
p.fNameTm3 = 'Tm3';
p.fNameCT1 = 'CT1';

% Conductance gains
p.gMi9 = 0.3;
p.gMi1 = 0.1 / 2;
p.gTm3 = 0.1 / 2;
p.gMi4 = 0.3 / 2;
p.gCT1 = 0.3 / 2;

%% Simulation parameters

% Number of realizations
p.numRep = 1000;

% Number of bootstraps
p.nboot = 1000;

% Batch size
p.batchSize = 100;

% Kernel extraction method
p.olsKernel = false;

%% Natural scene simulation parameters

% Standard deviation of Gaussian distribution of velocities
p.velStd = 100;

%% Extract user inputs

% Check for input arguments
if nargin > 0
    
    % Validate the number of input arguments
    if rem(nargin,2) ~= 0
        error('Arguments must be in name-value pairs!');
    end
    
    % Extract name-value pairs and store in parameter structure
    for ind = 1:2:nargin
        p.(varargin{ind}) = varargin{ind+1};
    end
end

%% Make temporal + spatial vectors and masks for stimulus generation

% Total simulation time
p.tTot = 2*p.tInt + p.tOn;

% Compute spatial position vector
p.xTot = p.xExtent;
p.x = (0:p.dx:p.xExtent-p.dx);

% Compute temporal position vector
if p.tInt <= 0
    p.t = (0:p.dt:p.tOn-p.dt)';
else
    p.t = (-p.tInt:p.dt:p.tOn+p.tInt-p.dt)';
end

% Define stimulus presentation mask
if p.tInt <= 0
    p.mask = ones(size(p.t), 'logical');
else
    p.mask = (p.t>=0) & (p.t<p.tOn);
end

% Define averaging mask
if p.tAv <= 0
    p.averagingMask = p.mask;
else
    p.averagingMask = (p.t >= p.tAv) & (p.t<p.tOn);
end

end

