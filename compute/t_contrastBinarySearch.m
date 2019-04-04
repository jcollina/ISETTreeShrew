%% t_contrastBinarySearch
% Use a binary SVM and binary search algorithm to create a contrast
% sensitivity function for the tree shrew.
%
% Description:
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
compute = 0;

loadData = 1;

toDivide = 20;

psfSigma = 12;

%note:usually 10/100 seconds, experimenting here
integrationTime = 50/1000;

sizeDegs = 5;

nTrialsNum = 1000;

% If you're not interested in computing a new CSF, but instead want to plot
% prevously gathered data, then set compute = 0 and specify the data you
% want to load. The data needs to have the following
oldDataToLoad = 'csf_1000_trials_psf_12_size_5.mat';

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
    contrastRange = [.001,.06];
    
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
            [acc,prec,t] = getSVMAcc(theMosaic, testScene, nullScene, nTrialsNum,psfSigma);
            
            % Save variables
            contrasts(cycle+1) = currentContrast;
            accuracies(cycle+1) = acc;
            
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
                finalAcc = acc;
                finalPrec = prec;
                %disp('next :)')
                break;
            end
            
            % If you don't make it, find the closest contrast and take that
            if cycle > maxCycles
                temp = [abs(75-rmmissing(accuracies)) ; rmmissing(accuracies) ; rmmissing(contrasts)]';
                temp = sortrows(temp,1);
                
                thresholdContrast = temp(1,3);
                finalAcc = temp(1,2);
                %disp('not quite there :/')
                break
                
            end
            
            cycle = cycle + 1;
            
        end
        
        % Save variables
        contrastsTotal{1,i} = rmmissing(contrasts);
        accuraciesTotal{1,i} = rmmissing(accuracies);
        thresholdContrasts(1,i) = thresholdContrast;
        finalAccuracy(1,i) = finalAcc;
    end
    time = toc(T);
    
    if ~exist('psfSigma','var')
        psfSigma = 7;
    end
    dataToSave = sprintf('csf_%.0f_trials_psf_%.0f_size_%.0f.mat',nTrialsNum,psfSigma,sizeDegs);
    if ~exist(dataToSave, 'file')
    save(strcat(dataToSave,'2'),'theMosaic','contrastsTotal','accuraciesTotal','thresholdContrasts','finalAccuracy','frequencyRange','nTrialsNum','time','sizeDegs','psfSigma','integrationTime') ;
    else
    end
    dataToLoad = dataToSave;
elseif loadData
    data = load(oldDataToLoad);
end

%% Plot data
%
% Visualize relationship between contrast and accuracy

plotBinarySearch(data)

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
