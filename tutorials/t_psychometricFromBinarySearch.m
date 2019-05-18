function t_psychometricFromBinarySearch(varargin)
%
%% t_psychometricFromBinarySearch (Human/Treeshrew)
% Use sampledThresholdBinarySearch, getPsychometricFit amd plotPsychometricFunction 
% to confirm that binay searches actually estimate the threshold

% Syntax:
%   t_psychometricFromBinarySearch('csfData','load('bestCSFData.mat'))
%   t_BinarySearchCSF('dataName','bestDataEver')
%   t_BinarySearchCSF('csfDataToPlot',load('treeshrew_csf_250_trials.mat'))
%

% Description:
% In Casagrande's 1984 treeshrew CSF paper, the researchers fit their data
% of contrast vs accuracy with a psychometric function, and used it to find
% their desired contrast threshold.
%
% In ISETBio, we don't have an a priori idea of the threshold location, so
% we set a wide contrast range to search. If we evenly sampled the full
% range, we would need a very large number of data points in order to fully
% sample the psychometric function. Instead, we use a binary search
% algorithm to find the desired contrast threshold, as can be seen in
% t_BinarySearchCSF.
%
% Nevertheless, we assume that there is an underlying psychometric
% function. And once we've found the contrast threshold, we can sample
% around that threshold and fit the local data to a psychometric function,
% proving that our binary search does, in fact, approximate the correct
% contrast values.

% Optional key/value pairs:
%   csfData - Binary search data (structure)
%   psychData - Previously gathered psychometric function data, if you just
%       want to plot (structure)
%   dataName - If you have a specific choice for the name of the generated
%       data, you can input dataName as a character array.
%   overwrite - If there is already a dataset with the same name, do you
%       want to overwrite it? true/false.
%   saveToWS - If you're just exploring the code, you may want to save the
%       generated data to the workspace after it is created. Logical.

% See also: t_BinarySearchCSF
%

% Tools used: sampledThresholdBinarySearch, getPsychometricFit,
% plotPsychometricFunction
%

% History:
%   04/02/19 jsc  Wrote initial version.

% Required csfData dataset can be created using t_BinarySearchCSF. Data must
% have the variables:
%       species
%       expInfo: struct w/
%           theMosaic
%           theOI
%           nTrialsNum
%           nFolds
%           psfSigma
%       binaryResults: struct w/
%           contrasts
%           accuracies
%           frequencyRange

% What data do you want to use? This could be overwritten by function
% input.

data = load('sampleCSFResults.mat');

% Parse optional input from function
p = inputParser;
p.StructExpand = false;
p.addParameter('csfData', data, @isstruct)
p.addParameter('psychData', struct, @isstruct)
p.addParameter('dataName', char.empty, @ischar)
p.addParameter('overwrite', false, @islogical)
p.addParameter('saveToWS', false, @islogical)


p.parse(varargin{:});
data = p.Results.csfData;
psychData = p.Results.psychData;
dataName = p.Results.dataName;
overwrite = p.Results.overwrite;
saveToWS = p.Results.saveToWS;

% Assume that we are creating a new dataset, unless indicated otherwise by
% passing a previously created psych dataset into the function.
compute = true;

if ~ismember('psychData',p.UsingDefaults)
    compute = false;
    if ~isempty(dataName)
        warning('You are loading data instead of generating new data, so no new data will be saved under your chosen name.')
    end
end

if ~ismember('csfData',p.UsingDefaults) && compute == false
    error('Either choose a psychometric function data structure to plot, or create one from binary search data.')
end

%% Parameters to determine which search to use
%
% In the binary search dataset, there will be multiple searches, for
% different
% spatial frequencies. We want to find a search that fits a few criteria in
% order to ensure we have enough different data points, and that the search
% matches set standards.

%First, we don't want to sample the entire field, and want to take advantage of
%the binary search we've already done. How many steps do we want to take
%before deciding we're close enougn to the threshold?

% 1 : 1/2 of the original contrast range
% 2 : 1/4 of the original contrast range
% ...
stepsBeforePlotting = 3;

% What's the minimum number of unique points you want to sample? So, the
% points after the previous parameter.
minInputPoints = 2;

% What's the range of satisfactory accuracies you want to have reached in the original
% search?
acceptedAccRange = [74.5,75.5];

% How many points within the threshold range do you want to
% compute?
numOutputPoints = 5;

% How many points outside the threshold range would you like to compute?
% These will be evenly spaced, spanning the range of contrasts searched.
numOtherPoints = 5;

if compute
    
    thresholdSample = sampledThresholdBinarySearch(data, ...
        'stepsBeforePlotting',stepsBeforePlotting, ...
        'acceptedAccRange',acceptedAccRange, ...
        'minInputPoints',minInputPoints, ...
        'numOtherPoints',numOtherPoints, ...
        'numOutputPoints',numOutputPoints ...
        );
    
    %% Fit the data with a psychometric function
    %
    % Determine number of data points used to compute each accuracy
    % measure
    expInfo = data.expInfo;
    nTrialsNum = expInfo.nTrialsNum;
    nTrialsPerMeasure = nTrialsNum/expInfo.nFolds;
    % Fit the data
    psychFitResults = getPsychometricFit(thresholdSample.contrasts,thresholdSample.accuracies, nTrialsPerMeasure ,'accThreshold',mean(acceptedAccRange));
    
    %% Save the data, including the results of the fit
    
    % Will data be saved to the workspace?
    if saveToWS
        assignin('base','thresholdSample',thresholdSample)
        assignin('base','psychFitResults',psychFitResults)
    end
    
    % Will the data be saved under a generic or chosen name?
    if isempty(dataName)
        psychDataToSave = sprintf('psychFn_%.0f_trials.mat',nTrialsNum);
    elseif endsWith(dataName,'.mat')
        psychDataToSave = dataName;
    else
        psychDataToSave = strcat(dataName,'.mat');
    end
    
    % Will the data be overwritten?
    if ~overwrite
        k = 1;
        psychDataToSaveTemp = psychDataToSave;
        while exist(psychDataToSaveTemp, 'file')
            psychDataToSaveTemp = strcat('v',num2str(k),'_',psychDataToSave);
            k = k + 1;
        end
        psychDataToSave = psychDataToSaveTemp;
    end
    
    save(psychDataToSave,'expInfo','thresholdSample','psychFitResults')
else
    % If you passed in a previously gathered dataset with the proper fields,
    % we'll expand that here.
    thresholdSample = psychData.thresholdSample;
    psychFitResults = psychData.psychFitResults;
end

%% Plot the psychometric function and the contrast threshold
%
% Either way, plot the data
plotPsychometricFunction(thresholdSample, 'fitResults', psychFitResults)

end