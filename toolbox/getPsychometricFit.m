function fitResults = getPsychometricFit(contrasts,percentCorrect,nTrials,varargin)

p = inputParser;

p.addParameter('accThreshold', 75, @isnumeric)

p.parse(varargin{:});

accThreshold = p.Results.accThreshold;

% Set up psychometric function model. Here we use a cumulative Weibull function
psychometricFunctionModel = @PAL_Weibull;

fractionCorrect = percentCorrect/100;

% Set up search grid
gridLevels = 100;
searchGridParams.alpha = logspace(log10(min(contrasts)),log10(max(contrasts)),gridLevels);
searchGridParams.beta = 10.^linspace(-4,4,gridLevels);
searchGridParams.gamma = 0.5;
searchGridParams.lambda = 0.0;

% Optimization settings for the fit
optionsParams             = optimset('fminsearch');
optionsParams.TolFun      = 1e-09;
optionsParams.MaxFunEvals = 1000;
optionsParams.MaxIter     = 1000;
optionsParams.Display     = 'off';

% Parameters for the curve fitting
% Parameters that are allowed to vary
% The parameters are: threshold, slope, guess-rate, lapse-rate
paramsFree = [1 1 0 0];
trialsNumCorrectPerContrastLevel = round(nTrials*fractionCorrect);
trialsNumPerContrastLevel = repmat(nTrials,1,length(fractionCorrect));

%% Fit the data and get the best fit params
paramsValues = PAL_PFML_Fit(contrasts(:), trialsNumCorrectPerContrastLevel(:), trialsNumPerContrastLevel(:), ...
            searchGridParams, paramsFree, psychometricFunctionModel, 'SearchOptions', optionsParams);

% Obtain the threshold at which performance cross a threshold performance, here 75%
performanceThreshold = accThreshold;
fitResults.contrastThreshold = psychometricFunctionModel(paramsValues, performanceThreshold, 'inverse');

%
% Obtain a high resolution version of the fitted function
fitResults.hiResContrasts = searchGridParams.alpha;
fitResults.hiResPerformance = 100 * PAL_Weibull(paramsValues, fitResults.hiResContrasts);

end