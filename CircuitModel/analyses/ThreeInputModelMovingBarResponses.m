function [ params, filters, meanResp, meanNumResp, meanDenResp, meanNumLnResp ] = ThreeInputModelMovingBarResponses(v, config, barParam, varargin)

%% Set parameters

% Set overall parameters
[ params ] = SetModelParameters(varargin{:});

% Make the filters
[ filters ] = MakeModelFilters(config, params);

% Define stimulus names
legendStr = {'PD ON','ND ON','PD OFF','ND OFF'};

%% Compute responses

% Compute responses
tic;
numV = length(v);
meanResp = nan(numV,4);

if nargout > 5
    meanNumResp = nan(numV,4);
    meanDenResp = nan(numV,4);
    meanNumLnResp = nan(numV, 4);
end

% Iterate over velocities
tAll = tic;
for indV = 1:numV
    
    % Make the stimulus array
    [ stimArray ] = MovingBars(params, barParam, v(indV));
    %     [ stimArray ] = MovingBarsWithStaticBarInterleave(params, barParam, v(indV));
    
    % Determine whether component responses need to be computed...
    if nargout > 5
        if params.useMeasuredFilters
            [ meanResp(indV,:), ~, ~, ...
                meanNumResp(indV,:), ...
                meanDenResp(indV,:),...
                meanNumLnResp(indV,:)...
                ] = ComputeSynapticModelResponseWithMeasuredFilters(stimArray, params, filters);
            
        else
            [ meanResp(indV,:), ~, ~, ...
                meanNumResp(indV,:), ...
                meanDenResp(indV,:),...
                meanNumLnResp(indV,:)...
                ] = ComputeSynapticModelResponse(stimArray, params, filters);
        end
    else
        if params.useMeasuredFilters
            [ meanResp(indV,:) ] = ComputeSynapticModelResponseWithMeasuredFilters(stimArray, params, filters);
        else
            [ meanResp(indV,:) ] = ComputeSynapticModelResponse(stimArray, params, filters);
        end
    end
    
end

fprintf('Completed simulation in %f seconds.\n', toc(tAll));

end

