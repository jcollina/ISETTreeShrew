function [percentCorrect,time] = getSVMAcc(theMosaic, testScene, nullScene, varargin)
% Calculate the accuracy of an SVM trained to discriminate between a test
% stimulus and a null stimulus using a given cone mosaic.
%
% Syntax:
%   [percentCorrect,time] = getSVMAcc(theMosaic, testScene, nullScene, varargin)
%

% Optional arguments:
% nTrialsNum:   Number of iterations to be trained upon
% species:      'treeshrew' or 'human' optics?
% nFolds:       How many folds to be used in k-fold validation?
% psfSigma:     For the optics, what's the spread of the point-spread
%               function?

p = inputParser;
p.addParameter('nTrialsNum', 250, @isnumeric);
p.addParameter('species', 'treeshrew', @ischar);
p.addParameter('nFolds', 10, @isnumeric);
p.addParameter('psfSigma', 7, @isnumeric);

% Parse input
p.parse(varargin{:});

psfSigma = p.Results.psfSigma;
nFolds = p.Results.nFolds;
nTrialsNum = p.Results.nTrialsNum;
species = p.Results.species;

tic;

emPathLength = 1;
taskIntervals = 2;
emPath = zeros(nTrialsNum, emPathLength, 2);

% Generate wavefront-aberration derived ts optics
switch species
    case 'treeshrew'        
        theOI = oiTreeShrewCreate(...
            'inFocusPSFsigmaMicrons', psfSigma ...
            );
    case 'human'
        theOI = oiCreate(...
            'inFocusPSFsigmaMicrons', psfSigma ...
            );
    otherwise
        error('species should be treeshrew or human')
end
% Compute the retinal image of the test stimulus
theTestOI = oiCompute(theOI, testScene);

% Compute the retinal image of the null stimulus
theNullOI = oiCompute(theOI, nullScene);

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
[pcVectors, ~, ~, ~, ~] = pca(classificationMatrix);

% Project the responses onto the space formed by the first 2 PC vectors
pcComponentsNum = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNum);

% Train a binary SVM classifier
svm = fitcsvm(classificationMatrixProjection,classes);

% Perform a 10-fold cross-validation on the trained SVM model
CVSVM = crossval(svm,'KFold',nFolds);

% Compute classification loss for the in-sample responses using a model
% trained on out-of-sample responses
%
fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
%
% Percent correct across all folds
percentCorrect = fractionCorrect.*100;%mean(fractionCorrect) * 100;

time = toc;

end
