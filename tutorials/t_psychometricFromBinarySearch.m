function t_psychometricFromBinarySearch(varargin)
%% t_psychometricFromBinarySearch (Human/Treeshrew)
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

% See also: t_BinarySearchCSF
%

% Tools used: getSVMAcc, getPsychometricFit
%

% History:
%   04/02/19 jsc  Wrote initial version.

% Load dataset created using t_BinarySearchCSF. Data must have the
% variables:
%       theMosaic
%       theOI
%       nTrialsNum
%       psfSigma
%       binaryResults: struct w/
%           contrastsTotal
%           frequencyRange

% What data do you want to use? This would be overwritten by function
% input.

data = load('sampleCSFData.mat');

% Parse optional input from function
p = inputParser;
p.StructExpand = false;
p.addParameter('csfData', data, @isstruct)
p.addParameter('psychData', struct, @isstruct)
p.addParameter('dataName', char.empty, @ischar)
p.addParameter('overwrite', false, @islogical)

p.parse(varargin{:});
data = p.Results.csfData;
psychData = p.Results.psychData;
dataName = p.Results.dataName;
overwrite = p.Results.overwrite;

% Load a previously generated set of data, if indicated.

% Assume that we are creating a new dataset, unless indicated otherwise by
% passing a previously created psych dataset into the function.
compute = true;

% Load a previously generated set of data, if indicated.
if ~ismember('psychData',p.UsingDefaults)
    compute = false;
    if ~isempty(dataName)
        warning('You are loading data instead of generating new data, so no new data will be saved under your chosen name.')
    end
end

if ~ismember('csfData',p.UsingDefaults) && compute == false
    error('Either choose a psychometric function data structure to plot, or create one from binary search data.')
end


nTrialsNum = data.expInfo.nTrialsNum;
nFolds = 10;%data.expInfo.nFolds;

%sizeDegs = round(theMosaic.fov);

% All we need now is the results of the binary search
%data = data.binaryResults;

%}

%% Selecting spatial frequency to use
%
% In the binary search dataset, you will have run the search for multiple
% spatial frequencies. We want to find a search that fits two criteria:
%   1)  We want the search to have a minimum number of steps (minSteps), to
%       ensure that there are multiple points close to the threshold.
%   2)  In the last four contrasts sampled, we want there to be at least
%       three unique values, so that we get an idea of the slope of the
%       psychometric function, not just a single point.

%We don't want to sample the entire field, and want to take advantage of
%the binary search we've already done. How many steps do we want to take
%before deciding we're close enougn to the threshold?

% 1 : 1/2 of the original contrast range
% 2 : 1/4 of the original contrast range
% ...
stepsBeforePlotting = 3;

% What's the minimum number of unique points you want to sample?
minInputPoints = 5;

% What's the range of satisfactory accuracies you want to have reached in the original
% search?
accRange = [74.5,75.5];

% How many points within the threshold range do you want to
% compute?
numOutputPoints = 5;

% How many points outside the threshold range would you like to compute?
% These will be evenly spaced, spanning the range of contrasts searched.
numOtherPoints = 5;

if compute

    sampledData = sampleThresholdBinarySearch(data, ...
        'stepsBeforePlotting',stepsBeforePlotting, ...
        'accRange',accRange, ...
        'minInputPoints',minInputPoints, ...
        'numOtherPoints',numOtherPoints, ...
        'numOutputPoints',numOutputPoints ...
        );

    % Fit the data with a psychometric function
    psychFitResults = getPsychometricFit(sampledData.contrasts,sampledData.accuracies, nTrialsNum/nFolds ,'accThreshold',mean(accRange));
    % contrasts, accuracies, threshold, ntrials
    % Save the data, including the results of the fit
    
    if isempty(dataName)
        psychDataToSave = [species,sprintf('psychFn_%.0f_trials.mat',nTrialsNum)];
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
    
    thresholdSample.frequency = spatFreq;
    thresholdSample.contrastsTotal = contrastsToFit;
    thresholdSample.accuraciesTotal = accuraciesToFit;
    thresholdSample.totalSE = acc_SE;
    
    save(psychDataToSave,'nTrialsNum','thresholdSample','psychFitResults')
else
    thresholdSample = psychData.thresholdSample;
    psychFitResults = psychData.psychFitResults;
end

%% Psychometric function fit implementation
%
% Plot the simulated data, the psychometric function, and the reported contrast threshold
plotPsychometricFunction(thresholdSample,psychFitResults)

%% Functions
    

end