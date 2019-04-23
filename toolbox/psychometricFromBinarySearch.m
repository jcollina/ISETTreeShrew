function [contrastsToPlot,accuraciesToPlot,contrastThreshold,hiResContrasts,hiResPerformance] = psychometricFromBinarySearch(data)

sizeDegs = data.sizeDegs;
theMosaic = data.theMosaic;
nTrialsNum = data.nTrialsNum;
psfSigma = data.psfSigma;

%
% %at this point, we have a psychometric function of accuracy as a function
% %of contrast that should sample the threshold well

% if I want to, I can compute an actual psychometric function here:

%pick the first one I chose
%%
j = 1;

numSearchPoints = 10;

nTrials = nTrialsNum/10;

while 1
    
    a = data.contrastsTotal{1,j};
    b = a(length(a)-numSearchPoints:length(a));
    
    if a > 12 && length(unique(b)) > 5
        break;
    end
    j = j+1;
end

contrastsTemp = data.contrastsTotal{1,j};

spatFreq = data.frequencyRange(j);

thresholdContrasts = contrastsTemp(length(contrastsTemp) - ...
    numSearchPoints:length(contrastsTemp));

minC = min(thresholdContrasts);
maxC = max(thresholdContrasts);

thresholdSamples = minC:(maxC-minC)/5:maxC;

contrastsToPlot = sort([.001,.005,.01,.015,.02,.03,thresholdSamples]);
%%
%%

% Size of testing set, based on 10-fold cross-evaluation

% Create presentation display
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);

% Create parameter structure for a low spatial frequency Gabor stimulus
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', data.frequencyRange(j), ...  % changing cycles/deg
    'orientationDegs', 0, ...               % 0 degrees
    'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
    'sizeDegs', data.sizeDegs, ...               % 14 x 14 size
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

parfor i = 1:length(contrastsToPlot)
    
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', spatFreq, ...  % changing cycles/deg
        'orientationDegs', 0, ...               % 0 degrees
        'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
        'sizeDegs', sizeDegs, ...               % 14 x 14 size
        'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
        'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
        'pixelsAlongWidthDim', [], ...          % pixels- width dimension
        'pixelsAlongHeightDim', [] ...          % pixel- height dimension
        );
    
    stimParams.contrast = contrastsToPlot(i);
    
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);
    
    [acc,~,~] = getSVMAcc(theMosaic, testScene, nullScene, nTrialsNum,psfSigma);
    
    accuraciesToPlot(i) = acc;
    
end
toc
%%
[contrastThreshold,hiResContrasts,hiResPerformance] = getPsychometricFit(contrastsToPlot,accuraciesToPlot/100,nTrials);
%%
figure()
plot(hiResContrasts,100*hiResPerformance, 'r-', 'LineWidth', 1.5);

hold on

plot(contrastsToPlot,accuraciesToPlot, 'ko', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);

line([contrastThreshold,contrastThreshold],[40,75])
line([min(contrastsToPlot),contrastThreshold],[75,75])
%
if max(contrastsToPlot) > .031
    set(gca,'xlim',[min(contrastsToPlot),max(contrastsToPlot)])
    contrastTicks = [0.005,0.01 0.02 0.03, max(contrastsToPlot) ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', int2str(round(maxCont))};
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
