function t_psychometricFromBinarySearch
%% t_psychometricFromBinarySearch (Human/Treeshrew)
% In Casagrande's 1984 treeshrew CSF paper, the researchers fit their data
% of contrast vs accuracy with a psychometric function, and used it to find
% their desired contrast threshold.
%
% In ISETBio, we don't have an a priori idea of the threshold location, so
% we set a wide contrast range to search. If we plotted evenly spaced
% points over this range, we would be unlike to sample the accuracies
% between ceiling and randomness, at least not enough to provide evidence
% for a psychometric function fit. Instead, we use a binary search
% algorithm to find the desired contrast threshold, as can be seen in
% t_BinarySearchCSF.
%
% Nevertheless, we assume that there is an underlying psychometric
% function. And once we've found the contrast threshold, we can sample
% around that threshold and fit the local data to a psychometric function,
% proving that our binary search does, in fact, locate the correct contrast
% values.

% See also: t_BinarySearchCSF
%

% Tools used: getSVMAcc, getPsychometricFit
%

% History:
%   04/02/19 jsc  Wrote initial version.

% Load dataset created using t_BinarySearchCSF. Data must have the
% variables:
%       sizeDegs
%       theMosaic
%       nTrialsNum
%       psfSigma
%       contrastsTotal
%       frequencyRange

% Load the dataset, choosing the name that matches your simulated data
data = load('2cd_max_csf_1000_trials_psf_12_size_5.mat');

sizeDegs = data.sizeDegs;
theMosaic = data.theMosaic;
nTrialsNum = data.nTrialsNum;
psfSigma = data.psfSigma;

%% Selecting spatial frequency to use
%
% In the binary search dataset, you will have run the search for multiple
% spatial frequencies. We want to find a search that fits two criteria:
%   1)  We want the search to have a minimum number of steps (minSteps), to
%       ensure that there are multiple points close to the threshold.
%   2)  In the last four contrasts sampled, we want there to be at least
%       three unique values, so that we get an idea of the slope of the
%       psychometric function, not just a single point.

minSteps = 8;

j = 1;
while 1
    a = data.contrastsTotal{1,j};
    if length(a) > minSteps
        thresholdContrastsTemp = round(a(length(a)-4:length(a)),4);
        if length(unique(thresholdContrastsTemp)) > 3
           thresholdContrasts = unique(thresholdContrastsTemp);
           break;
        end
    end
    j = j+1;    
end

spatFreq = data.frequencyRange(j);

%contrastsTemp(length(contrastsTemp) - ...
%    numSearchPoints:length(contrastsTemp));

minC = min(thresholdContrasts);
maxC = max(thresholdContrasts);

numSamplePoints = 5;

thresholdSamples = minC:(maxC-minC)/numSamplePoints:maxC;

contrastsToPlot = sort([.001,.005,.01,.015,.02,.03,thresholdSamples]);

%% Using SVM to calculate accuracy for each contrast

% Size of testing set, based on 10-fold cross-evaluation

% Create presentation display
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);

% Create parameter structure for a low spatial frequency Gabor stimulus
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', data.frequencyRange(j), ...  % changing cycles/deg
    'orientationDegs', 0, ...               % 0 degrees
    'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
    'sizeDegs', data.sizeDegs, ...          % size
    'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
    'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
    'pixelsAlongWidthDim', [], ...          % pixels- width dimension
    'pixelsAlongHeightDim', [] ...          % pixel- height dimension
    );

% Generate null scene (contrast = 0)
stimParams.contrast = 0.0;

nullScene = generateGaborScene(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay);

accuraciesToPlot = zeros(1,length(contrastsToPlot));

% Cycle through contrasts, using the SVM to calculate accuracy for each. 

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,length(contrastsToPlot)) '\n\n']);

parfor i = 1:length(contrastsToPlot)
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', spatFreq, ...  % changing cycles/deg
        'orientationDegs', 0, ...               % 0 degrees
        'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
        'sizeDegs', sizeDegs, ...               % size
        'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
        'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
        'pixelsAlongWidthDim', [], ...          % pixels- width dimension
        'pixelsAlongHeightDim', [] ...          % pixel- height dimension
        );
    
    % Set the contrast level for this iteration
    stimParams.contrast = contrastsToPlot(i);
    
    % Generate a Gabor patch with that contrast
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay,...
        'minimumPixelsPerHalfPeriod', 5);
    
    % Determine how well the specified optics can discriminate between that
    % patch and a null stimulus
    [acc,~] = getSVMAcc(theMosaic, testScene, nullScene, 'nTrialsNum',nTrialsNum,'psfSigma',psfSigma);%,'species',species);
    accuraciesToPlot(i) = mean(acc);
    fprintf('\b|\n');    
end
%% Psychometric function fit implementation
%
% Fit the data with a psychometric function
[contrastThreshold,hiResContrasts,hiResPerformance] = getPsychometricFit(contrastsToPlot,accuraciesToPlot/100,nTrialsNum/10);

% Then, plot the simulated data, the psychometric function, and the
% reported contrast threshold
plotPsychometricFunction(hiResContrasts,hiResPerformance,contrastsToPlot,accuraciesToPlot,contrastThreshold,spatFreq)

%% Functions
function plotPsychometricFunction(hiResContrasts,hiResPerformance,contrastsToPlot,accuraciesToPlot,contrastThreshold,spatFreq)
figure()
plot(hiResContrasts,100*hiResPerformance, 'r-', 'LineWidth', 1.5);
hold on
plot(contrastsToPlot,accuraciesToPlot, 'ko', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
line([contrastThreshold,contrastThreshold],[40,75])
line([min(contrastsToPlot),contrastThreshold],[75,75])
if max(contrastsToPlot) > .031
    set(gca,'xlim',[min(contrastsToPlot),max(contrastsToPlot)])
    contrastTicks = [0.005,0.01 0.02 0.03, max(contrastsToPlot) ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', int2str(round(maxContrastsToPlot))};
else
    set(gca,'xlim',[min(contrastsToPlot),.03])
    contrastTicks = [0.005 0.01 0.02 0.03];
    contrastTickLabels = {'0.005', '.01',  '.02',  '.03'};
end
set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);
set(gca, 'YLim', [40 105], 'XScale', 'log')
set(gca, 'FontSize', 16)
xlabel('\it Contrast (Michelson)');
ylabel('\it SVM Accuracy');
title(sprintf('Individual Psychometric Function for\nSpatial Frequency of %.0f cpd',spatFreq))
hold off
end

end