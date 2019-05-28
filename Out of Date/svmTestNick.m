
%simple code to check something that has been bugging me: does nTrials
%affect the accuracy past what we'd expect in terms of coarse/fine?

clear; close all;
%ieInit;

% number of iterations of the code, for looking at std error 
numIter = 2;
numTrials = 2;

percentCorrect = zeros(numIter,numTrials);

emPathLength = 1;

taskIntervals = 2;

theMosaic = coneMosaicTreeShrewCreate(75, ...%theOI.optics.micronsPerDegree, ...
    'fovDegs', 1, ...        % match mosaic width to stimulus size
    'integrationTimeSeconds', 10/1000);

%% Create presentation display and place it 5 cm in front of the eye
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);


%picked a contrast that has given non ceiling/floor accuracies in the past
contrast = .003;

% parameter struct for a low spatial frequency Gabor stimulus
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 1, ... % changing cycles/deg
    'orientationDegs', 0, ...               % 0 degrees
    'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
    'sizeDegs', 1, ...                     % 14 x 14 size
    'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
    'contrast', contrast,...                   % 0.005 Michelson contrast
    'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
    'pixelsAlongWidthDim', [], ...          % pixels- width dimension
    'pixelsAlongHeightDim', [] ...          % pixel- height dimension
    );

%
% Generate a scene representing the Gabor stimulus with the above params as
% realized on the presentationDisplay
testScene = generateGaborScene(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay);

% Visualize different aspects of the generated scene
visualizeScene(testScene, 'displayRadianceMaps', false);
%
% zero contrast for the null stimulus
stimParams.contrast = 0.0;

nullScene = generateGaborScene(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay);




% Generate wavefront-aberration derived ts optics
theOI = oiTreeShrewCreate();

% Compute the retinal image of the test stimulus
theTestOI = oiCompute(theOI, testScene);

% Compute the retinal image of the null stimulus
theNullOI = oiCompute(theOI, nullScene);

for i = 1:numIter
for j = 1:numTrials


nTrialsNum = 2^(2+j);
bigN(j) = nTrialsNum;

emPath = zeros(nTrialsNum, emPathLength, 2);

%
%
% Compute mosaic excitation responses to the test stimulus
coneExcitationsTest = theMosaic.compute(theTestOI, 'emPath', emPath);

% Compute mosaic excitation responses to the null stimulus
coneExcitationsNull = theMosaic.compute(theNullOI, 'emPath', emPath);

% Obtain the indices of the grid nodes that contain cones
[~,~,~, nonNullConeIndices] = theMosaic.indicesForCones;

% Extract the response vectors for nodes containing cones
[nTrials, nRows, mCols, nTimeBins] = size(coneExcitationsTest);
coneExcitationsTestReshaped = reshape(coneExcitationsTest, [nTrials nRows*mCols nTimeBins]);
coneExcitationsNullReshaped = reshape(coneExcitationsNull, [nTrials nRows*mCols nTimeBins]);
testResponses = coneExcitationsTestReshaped(:, nonNullConeIndices, :);
nullResponses = coneExcitationsNullReshaped(:, nonNullConeIndices, :);

% Collapse response vectors across space and time
responseSize = numel(nonNullConeIndices)*nTimeBins;
testResponses = reshape(testResponses, [nTrials responseSize]);
nullResponses = reshape(nullResponses, [nTrials responseSize]);

% Form the classification matrix.

%
if (taskIntervals == 1)
    % In the one interval task, the null and test response instances are labelled as the 2 classes.
    % Allocate matrices
    classificationMatrix = nan(2*nTrials, responseSize);
    classes = nan(2*nTrials, 1);
    % Class 1
    classificationMatrix(1:nTrials,:) = nullResponses;
    classes((1:nTrials)) = 0;
    % Class 2
    classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
    classes(nTrials+(1:nTrials)) = 1;
elseif (taskIntervals == 2)
    % In the two inteval task, we concatenate [null test] as one class and [test null] as the other.
    % Allocate matrices
    classificationMatrix = nan(nTrials, 2*responseSize);
    classes = nan(nTrials, 1);
    halfTrials = floor(nTrials/2);
    % Class 1
    classificationMatrix(1:halfTrials,:) = [...
        nullResponses(1:halfTrials,:) ...
        testResponses(1:halfTrials,:)];
    classes((1:halfTrials)) = 0;
    % Class 2
    idx = halfTrials+(1:halfTrials);
    classificationMatrix(idx,:) = [...
        testResponses(idx,:) ...
        nullResponses(idx,:)];
    classes(idx) = 1;
else
    error('Task can have 1 or 2 intervals only.')
end

% Find principal components of the responses
[pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

% Project the responses onto the space formed by the first 2 PC vectors
pcComponentsNum = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNum);

% Train a binary SVM classifier
svm = fitcsvm(classificationMatrixProjection,classes);

% Perform a 10-fold cross-validation on the trained SVM model
kFold = 10;
CVSVM = crossval(svm,'KFold',kFold);

% Compute classification loss for the in-sample responses using a model
% trained on out-of-sample responses
%
fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
%
% Average percent correct across all folds
percentCorrect(i,j) = mean(fractionCorrect) * 100;

end
end

sem = std(percentCorrect,0,1)./bigN;
%%
figure()
errorbar(bigN,mean(percentCorrect,1),sem, 's-')
set(gca, 'XLim', [0.5*bigN(1) 1.1*bigN(end)], 'YLim', [0 100], 'YTick', 0:20:100, 'XTick', [0:100:bigN(end)]);
grid on;
title('Effect of Changing nTrialsNum on SVM Accuracy')
xlabel('nTrialsNum')
ylabel('% Correct')
