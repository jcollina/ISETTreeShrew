function sampledData = sampledThresholdBinarySearch(samples,performance,,varargin)
% Sample the psychometric threshold space indicated by a binary search
%
% Syntax:
%     thresholdSample = sampledThresholdBinarySearch(data, ...
%         'stepsBeforePlotting',3, ...
%         'acceptedAccRange',[74.5,75.5], ...
%         'minInputPoints',5, ...
%         'numOtherPoints',5, ...
%         'numOutputPoints',5, ...
%         'contrastRange',[.001,.03] ...
%         );
%
% Inputs:
%   thresholdSample - Matlab structure detailing SVM accuracy as a function
%   of a stimulus feature. With fields

%   'contrasts', 'accuracies', 'frequency' and 'contrastRange'.
  
%
% Optional key/value pairs:
%   stepsBeforePlotting - steps to ignore in search in order to sample
%   threshold
%   acceptedAccRange - Not all searches reached the same accuracy. What
%   would you consider acceptable?
%   minInputPoints - Number of unique points after the stepsBeforePlotting
%   numOtherPoints - Number of points across the full contrast range to
%   also plot
%   numOutputPoints - Number of points to plot in the threshold

% See also: t_BinarySearchCSF, t_psychometricFromBinarySearch, plotPsychometricFunction

% History:
%   05/15/19 jsc  Wrote initial version.

% Parse inputs
p = inputParser;
p.addParameter('stepsBeforePlotting', 3, @isnumeric)
p.addParameter('minInputPoints', 5, @isnumeric)
p.addParameter('numOutputPoints', 5, @isnumeric)
p.addParameter('numOtherPoints', 5, @isnumeric)
p.addParameter('acceptedAccRange', [74.5 , 75.5],@isnumeric)
p.parse(varargin{:});
stepsBeforePlotting = p.Results.stepsBeforePlotting;
minInputPoints = p.Results.minInputPoints;
numOutputPoints = p.Results.numOutputPoints;
numOtherPoints = p.Results.numOtherPoints;
acceptedAccRange = p.Results.acceptedAccRange;


% Pull information about the simulated experiment to ensure the same
% parameters are used here
theMosaic = data.expInfo.theMosaic;
theOI = data.expInfo.theOI;
presentationDisplay = data.expInfo.presentationDisplay;
stimParams = data.expInfo.stimParams;

nTrialsNum = data.expInfo.nTrialsNum;
nFolds = data.expInfo.nFolds;


% Now, we need the results of the binary search
data = data.binaryResults;

% Search through the binary searches in order to find one that matches the
% requirements
j = 1;
while 1
    a = data.contrasts{1,j};
    thresholdContrastsTemp = unique(round(a(stepsBeforePlotting : length(a)),4));
    if acceptedAccRange(1) < data.finalAccuracy(j) && acceptedAccRange(2) > data.finalAccuracy(j)
        if length(thresholdContrastsTemp) > minInputPoints
            thresholdContrasts = thresholdContrastsTemp;
            break;
        end
    end
    j = j+1;
    if j > length(data.contrasts)
        error('Sadly, none of the binary searches in this dataset fulfill your requirements.')
    end
end

% Choose that frequency for the Gabor stimulus
stimParams.spatialFrequencyCyclesPerDeg = data.frequencyRange(j);

% Determine the contrasts to analyze
minC = min(thresholdContrasts)-range(thresholdContrasts)/2;
maxC = max(thresholdContrasts)+range(thresholdContrasts)/2;

thresholdSamples = minC:(maxC-minC)/numOutputPoints:maxC;
sampleContrastRange = min(data.contrastRange):range(data.contrastRange)/numOtherPoints:max(data.contrastRange);

contrasts = sort([sampleContrastRange,thresholdSamples]);

%% Using SVM to calculate accuracy for each contrast

% Size of testing set, based on 10-fold cross-evaluation


% Generate null scene (contrast = 0)
stimParams.contrast = 0.0;

nullScene = generateGaborScene(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay);

accuracies = zeros(1,length(contrasts));

%% Cycle through contrasts, using the SVM to calculate accuracy for each.

T = tic;

count = 1;
N = length(contrasts);

D = parallel.pool.DataQueue;

w = waitbar(0, 'Progress');
afterEach(D, @nUpdateWaitbar);

parfor i = 1:length(contrasts)
    
    % Set the contrast level for this iteration
    stimParamsTemp = stimParams;
    stimParamsTemp.contrast = contrasts(i);
    
    % Generate a Gabor patch with that contrastm
    testScene = generateGaborScene(...
        'stimParams', stimParamsTemp,...
        'presentationDisplay', presentationDisplay,...
        'minimumPixelsPerHalfPeriod', 5);
    
    % Determine how well the specified optics can discriminate between that
    % patch and a null stimulus
    acc = getSVMAcc(theMosaic, theOI, testScene, nullScene, 'nTrialsNum',nTrialsNum,'nFolds',nFolds);
    accuracies(i) = mean(acc);
    accuracySE(i) = std(acc)/sqrt(length(acc));
    
    send(D,i)
end

delete(gcp)
close(w)

sampledData.time = toc(T)/60;
sampledData.frequency = stimParams.spatialFrequencyCyclesPerDeg;
sampledData.contrasts = contrasts;
sampledData.accuracies = accuracies;
sampledData.accuracySE = accuracySE;
sampledData.acceptedAccRange = acceptedAccRange;

%% Functions

    function nUpdateWaitbar(~)
        waitbar(count/N, w);
        count = count + 1;
    end

end