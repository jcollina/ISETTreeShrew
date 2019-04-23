%% t_BinarySearchCSF
% Use a binary SVM and binary search algorithm to create a contrast
% sensitivity function.
%
% Description: 
%   Essentially trying to simulate Casagrande's 1984 CSF method in ISETBio.
%   The key steps in that paper are as follows:
%   1) Choose the specific spatial frequencies for which you will
%  determine the sensitivity to contrast
%   2) For each spatial frequency, determine the contrast level (aka 
% the contrast threshold) that causes the tree shrews to accurately 
% discern the contrast 75% of the time
%
%       ** In ISETBio, this step is accomplished by using a support vector
%       machine (SVM) approach to determine the fraction of the time that
%       the algorithm can detect a difference in the mosaic excitation
%       between a null scene and a scene with a contrast. A binary search
%       algorithm is used to identify the specific contrast that leads to
%       75% of the SVM predictions being accurate. 

%       The binary search works as follows: The spatial frequency of the
%       gabor stimulus is set at the beginning of the search, and is held
%       constant throughout. A minimum and maximum contrast is chosen, with
%       the assumption that the 75% contrast is within the chosen range. 
%       The steps involved in the binary search process are well-described
%       in https://en.wikipedia.org/wiki/Binary_search_algorithm.  

%       Because there is some variability in the SVM results, it is
%       possible for the binary search algorithm to focus on the wrong
%       contrast, without reaching the desired accuracy. For 1000 trials
%       and 10 folds, the SVM results usually have a standard error between
%       1% and 2%. We recognize this problem by setting a maximum number of
%       iterations. If that number of iterations is reached without the
%       desired contrast being found, then the code checks for the contrast
%       that led to an accuracy closest to 75%. If 75% falls within a 99%
%       CI for that SVM, then we use that contrast as our "threshold
%       contrast" and move on to the next spatial frequency.

%   3) Calculate the contrast sensitivity at each spatial frequency by
% taking the reciprocal of the threshold contrast 
%   4) Plot the sensitivity as a function of the spatial frequency

% ...
% See also:
%

% History:
%   03/14/19 jsc  Wrote initial version.

%% Initialize workspace and close old figures
%clear; close all;
ieInit;

%% What do you want to do?
%

% If you're not interested in computing a new CSF, but instead want to plot
% prevously gathered data, then set compute = 0 and specify the data you
% want to load.

compute = 0;

loadData = 1;
oldDataToLoad = '2cd_min_csf_1000_trials_psf_12_size_5.mat';

nameData = 0;
%dataName = 'maxConeCSF';

toDivide = 20;

%min cone separation of 6 -> ~ 32,000 cones/mm^2
%min cone separation of 7.5 -> ~ 22,000 cones/mm^2
%min cone separation of 8.5 -> ~ 16,000 cones/mm^2

cone_spacing = 6; %microns

psfSigma = 12;

%note:usually 10/1000 seconds, experimenting here
integrationTime = 10/1000;

sizeDegs = 5;

nTrialsNum = 1000;

psychFn = 0;

if compute
    
    %% Parameters
    %
    % How much data do you want to use to train the SVM?
    
    
    % What discrete sptial frequencies do you want to find the sensitivities
    % for?
    frequencyRange = 0.75:0.25:2;
    
    % What range of contrasts do you want to explore for the spatial
    % frequencies? Make it wide enough to include accuracies other than 100%
    % and 50%!
    contrastRange = [.001,.03];
    
    % For the binary search, what range do you want it to reach before
    % stopping?
    desiredAccRange = [74.5,75.5];
    
    % For the binary search, how many iterations do you want before stopping?
    % If you don't want to ever "give up" (only a good idea with a high
    % nTrialNum), set maxCycles very large.
    maxCycles = 15;
    
    % How big do you want the stimulus, and therefore the mosaic, to be?
    % This will dramatically affect the time it takes to run the code.
    
    % Create cone mosaic of the appropriate size
    theMosaic = coneMosaicTreeShrewCreate(75, ...
        'fovDegs', sizeDegs, ...
        'customLambda', cone_spacing, ...
        'integrationTimeSeconds', integrationTime);
    
    %%
    %% Initialize Variables
    %
    % Initialize cell array to store psychometric function data
    contrastsTotal = cell(1,length(frequencyRange));
    accuraciesTotal = cell(1,length(frequencyRange));
    
    % Initialize vector to store threshold contrasts for each spatial frequency
    thresholdContrasts = zeros(1,length(frequencyRange));
    finalAccuracy = zeros(1,length(frequencyRange));
    finalSE = zeros(1,length(frequencyRange));
    %% Main loop
    %
    % Begin looping over the discrete spatial frequencies chosen
    T = tic;
    parfor i = 1:length(frequencyRange)
        
        % Size of testing set, based on 10-fold cross-evaluation
        nTrials = nTrialsNum/10;
        
        % Create presentation display
        presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
        frequencyRange(i);
        % Create parameter structure for a low spatial frequency Gabor stimulus
        stimParams = struct(...
            'spatialFrequencyCyclesPerDeg', frequencyRange(i), ...  % changing cycles/deg
            'orientationDegs', 0, ...               % 0 degrees
            'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
            'sizeDegs', sizeDegs, ...               % 14 x 14 size
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
        
        %% Secondary loop
        
        % initialize parameters for binary search
        maxContrast = max(contrastRange);
        minContrast = min(contrastRange);
        contrasts = NaN(1,maxCycles);
        accuracies = NaN(1,maxCycles);
        standardError = NaN(1,maxCycles);
        cycle = 0;
        
        % Start search
        while 1
            
            currentContrast = mean([maxContrast,minContrast]);
            
            % Generate test stimulus based on current contrast of interest
            
            stimParams.contrast = currentContrast;
            
            testScene = generateGaborScene(...
                'stimParams', stimParams,...
                'presentationDisplay', presentationDisplay, ...
                'minimumPixelsPerHalfPeriod', 5);
            
            % Use a binary SVM to get a measure of accuracy in terms of determining
            % which stimulus had gratings
            %%
            [svm_results, t] = getSVMAcc(theMosaic, testScene, nullScene, 'nTrialsNum', 100, 'psfSigma' , psfSigma)
            %%
            acc = mean(svm_results);
            acc_SE = std(svm_results)/sqrt(length(svm_results)); 
            %%
            % Save variables
            contrasts(cycle+1) = currentContrast;
            accuracies(cycle+1) = acc;
            standardError(cycle+1) = acc_SE;
            
            
            % Determine parameters for next cycle of search
            if acc < min(desiredAccRange)
                minContrast = currentContrast;
                %disp('accuracy was low')
            elseif acc > max(desiredAccRange)
                maxContrast = currentContrast;
                %disp('accuracy was high')
            else
                % if reached desired range
                thresholdContrast = currentContrast;
                finalAcctemp = acc;
                finalSE = acc_SE;
                %disp('next :)')
                break;
            end
            
            % If you don't make it, find the closest contrast- if that's within 1, take that
            if cycle > maxCycles
                temp = [abs(75-rmmissing(accuracies)) ; rmmissing(accuracies) ; rmmissing(contrasts) ; rmmissing(standardError)]';
                temp = sortrows(temp,1);
                
                thresholdContrast = temp(1,3);
                tempAcc = temp(1,2);
                tempSE = temp(1,4);
                
                %if abs(75-tempAcc) < 1
                if (tempAcc-2*tempSE) < 75 && (tempAcc+2*tempSE) > 75
                %disp('not quite there :/')
                    finalAcctemp = tempAcc
                    finalSEtemp = tempSE
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
        
        % Save variables
        contrastsTotal{1,i} = rmmissing(contrasts);
        accuraciesTotal{1,i} = rmmissing(accuracies);
        thresholdContrasts(1,i) = thresholdContrast;
        finalAccuracy(1,i) = finalAcctemp;
        finalSE(1,i) = finalSEtemp;
    end
    time = toc(T);
    
    if ~exist('psfSigma','var')
        psfSigma = 7;
    end
    %%
    if nameDame
        dataToSave = dataName;
    else
        dataToSave = sprintf('csf_%.0f_trials_psf_%.0f_size_%.0f.mat',nTrialsNum,psfSigma,sizeDegs);
    end
    %%
    i = 1;
    
    dataToSaveTemp = dataToSave;
    
    while exist(dataToSaveTemp, 'file')
        dataToSaveTemp = strcat(num2str(i),dataToSave);
        i = i + 1;
    end
    %%
    save(dataToSaveTemp,'theMosaic','contrastsTotal','accuraciesTotal','thresholdContrasts','finalAccuracy','frequencyRange','nTrialsNum','time','sizeDegs','psfSigma','integrationTime','cone_spacing')
    %%
    data = load(dataToSaveTemp);
elseif loadData
    data = load(oldDataToLoad);
end

%% Plot data
%
% Visualize relationship between contrast and accuracy

plotBinarySearch(data)
%%
plotCSF(data,toDivide)

%% Functions

function plotBinarySearch(data)

maxCont = max(cellfun(@(x) max(x),data.contrastsTotal));
minCont = min(cellfun(@(x) min(x),data.contrastsTotal));

figure()
hold on
for i = 1:length(data.frequencyRange)
    tempMat = [data.contrastsTotal{1,i};data.accuraciesTotal{1,i}];
    mat = sortrows(tempMat');
    color = rand(1,3);
    plot(mat(:,1),mat(:,2),'LineWidth',2,'color',color)
    plot(data.thresholdContrasts(i),data.finalAccuracy(i),'*','MarkerSize',10,'color',color)
end

if maxCont > .031
    set(gca,'xlim',[minCont,maxCont])
    contrastTicks = [0.005,0.01 0.02 0.03, maxCont ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', sprintf('%.0f',maxCont)};
else
    set(gca,'xlim',[minCont,.03])
    contrastTicks = [0.005 0.01 0.02 0.03];
    contrastTickLabels = {'0.005', '.01',  '.02',  '.03'};
end
set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);

set(gca, 'XScale', 'log')
set(gca, 'FontSize', 16)
xlabel('\it Contrast (Michelson)');
ylabel('\it SVM Accuracy');
title('Psychometric Function Approximations')
hold off
end
