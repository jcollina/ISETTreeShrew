function [searchResults] = binarySearchOverContrast(contrastRange,theMosaic,theOI,varargin)

[defaultStimParams,defaultPresentationDisplay] = csfStimParamsDefault();

p = inputParser;
p.StructExpand = false;
p.addParameter('stimParams', defaultStimParams, @isstruct)
p.addParameter('presentationDisplay', defaultPresentationDisplay, @isstruct)
p.addParameter('maxCycles', 8, @isnumeric)
p.addParameter('desiredAccRange', [74.5 , 75.5],@isnumeric)
p.addParameter('frequency', 0, @isnumeric)
p.addParameter('nTrialsNum', 500, @isnumeric)
p.addParameter('nFolds', 10, @isnumeric)

p.parse(varargin{:});
stimParams = p.Results.stimParams;
presentationDisplay = p.Results.presentationDisplay;
maxCycles = p.Results.maxCycles;
desiredAccRange = p.Results.desiredAccRange;
nTrialsNum = p.Results.nTrialsNum;
nFolds = p.Results.nFolds;

if ~isfield(stimParams,'spatialFrequencyCyclesPerDeg')
    try
        stimParams.spatialFrequencyCyclesPerDeg = p.Results.frequency;
    catch
        warning('Need to select a spatial frequency!')
    end
end

% Generate null scene (contrast = 0)
stimParams.contrast = 0;

nullScene = generateGaborScene(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay);

maxContrast = max(contrastRange);
minContrast = min(contrastRange);
contrasts = NaN(1,maxCycles);
accuracies = NaN(1,maxCycles);
standardError = NaN(1,maxCycles);
cycle = 0;

while 1
    
    % Generate test stimulus based on current contrast of interest
    currentContrast = mean([minContrast,maxContrast]);
    
    stimParams.contrast = currentContrast;
    
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', 5);
    
    % Use a binary SVM to get a measure of accuracy in terms of
    % determining which stimulus had gratings
    svm_results = getSVMAcc(theMosaic, theOI, testScene, nullScene, 'nTrialsNum', nTrialsNum,'nFolds',nFolds);

    acc = mean(svm_results);
    acc_SE = std(svm_results)/sqrt(length(svm_results));

    % Save variables
    contrasts(cycle+1) = currentContrast;
    accuracies(cycle+1) = acc;
    standardError(cycle+1) = acc_SE;
    
    % Determine parameters for next cycle of search
    if acc < min(desiredAccRange)
        minContrast = currentContrast;
    elseif acc > max(desiredAccRange)
        maxContrast = currentContrast;
    else
        % if reached desired range
        thresholdContrast = currentContrast;
        finalAcc = acc;
        finalSEtemp = acc_SE;
        break;
    end
    
    % If you don't make it, find the closest contrast- if 75%
    % is within a 99% CI of that SVM, accept that contrast
    if cycle > maxCycles
        temp = [abs(75-rmmissing(accuracies)) ; rmmissing(accuracies) ; rmmissing(contrasts) ; rmmissing(standardError)]';
        temp = sortrows(temp,1);
        
        tempAcc = temp(1,2);
        tempSE = temp(1,4);
        
        if (tempAcc-1.6*tempSE) < 75 && (tempAcc+1.6*tempSE) > 75
            finalAcc = tempAcc;
            finalSEtemp = tempSE;
            thresholdContrast = temp(1,3);
            break
        else
            contrasts = NaN(1,maxCycles);
            accuracies = NaN(1,maxCycles);
            standardError = NaN(1,maxCycles);
            cycle = -1;
        end
        
    end
    cycle = cycle + 1;
end

% Create structure for results
searchResults.contrasts = rmmissing(contrasts);
searchResults.accuracies = rmmissing(accuracies);
searchResults.finalAcc = finalAcc;
searchResults.finalSE = finalSEtemp;
searchResults.thresholdContrast = thresholdContrast;
end