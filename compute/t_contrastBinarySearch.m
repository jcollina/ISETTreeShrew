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
clear; close all;
ieInit;

%% What do you want to do?
%
compute = 1;

% If compute == 0, and you
dataToLoad = 'dataTotal1000_.75-1.25_psf_15.mat';

toDivide = 12;

psychFn = 0;

if compute
    
    %% Parameters
    %
    % How much data do you want to use to train the SVM?
    nTrialsNum = 1000;
    
    % What discrete sptial frequencies do you want to find the sensitivities
    % for?
    frequencyRange = 0.75:0.25:2;
    
    % What range of contrasts do you want to explore for the spatial
    % frequencies? Make it wide enough to include accuracies other than 100%
    % and 50%!
    contrastRange = [.001,.05];
    
    % For the binary search, what range do you want it to reach before
    % stopping?
    desiredAccRange = [74.5,75.5];
    
    % For the binary search, how many iterations do you want before stopping?
    % If you don't want to ever "give up" (only a good idea with a high
    % nTrialNum), set maxCycles very large.
    maxCycles = 20;
    
    % How big do you want the stimulus, and therefore the mosaic, to be?
    % This will dramatically affect the time it takes to run the code.
    sizeDegs = 3;

    % Create cone mosaic of the appropriate size
    theMosaic = coneMosaicTreeShrewCreate(75, ...
        'fovDegs', sizeDegs, ...
        'integrationTimeSeconds', 10/1000);
    
%%
    %% Initialize Variables
    %
    % Initialize cell array to store psychometric function data
    contrastsTotal = cell(1,length(frequencyRange));
    accuraciesTotal = cell(1,length(frequencyRange));
    
    % Initialize vector to store threshold contrasts for each spatial frequency
    thresholdContrasts = zeros(1,length(frequencyRange));
    
    %% Main loop
    %
    % Begin looping over the discrete spatial frequencies chosen
    T = tic;
    for i = 1:length(frequencyRange)
        
        % Size of testing set, based on 10-fold cross-evaluation
        nTrials = nTrialsNum/10;
        
        % Create presentation display
        presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
        frequencyRange(i)
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
        accuracies = NaN(1,MaxCycles);
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
            [acc,t] = getSVMAcc(theMosaic, testScene, nullScene, nTrialsNum);
            
            % Save variables
            contrasts(cycle+1) = [contrasts, currentContrast];
            accuracies(cycle+1) = [accuracies, acc];
            
            % Determine parameters for next cycle of search
            if acc < min(desiredAccRange)
                minContrast = currentContrast;
                disp('accuracy was low')
            elseif acc > max(desiredAccRange)
                maxContrast = currentContrast;
                disp('accuracy was high')
            else
                % if reached desired range
                thresholdContrast = currentContrast;
                finalAcc = acc;
                disp('next :)')
                break;
            end
            
            % Allow for exit due to reaching maximum number of cycles
            if cycle > maxCycles
                thresholdContrast = currentContrast;
                exhaustedFinalAcc = acc;
                disp('next :/')
                break;
            end
            
            cycle = cycle + 1;
        end
        
        % Save variables
        contrastsTotal{1,i} = rmmissing(contrasts);
        accuraciesTotal{1,i} = rmmissing(accuracies);
        thresholdContrasts(1,i) = thresholdContrast;
    end
    time = toc(T);
    save('compute/dataTotal','theMosaic','contrastsTotal','accuraciesTotal','thresholdContrasts','frequencyRange','nTrialsNum','time','sizeDegs')    
else    
    load(dataToLoad)
    frequencyRange = [.75 1 1.25]
end

%% Plot data
%
% Visualize relationship between contrast and accuracy
figure(1)
hold on
xlim([.001,.03])
for i = 1:length(frequencyRange)
    tempMat = [contrastsTotal{1,i};accuraciesTotal{1,i}];
    mat = sortrows(tempMat');
    color = rand(1,3);
    plot(mat(:,1),mat(:,2),'LineWidth',2,'color',color)
    plot(tempMat(1,size(tempMat,2)),tempMat(2,size(tempMat,2)),'*','MarkerSize',10,'color',color)
end
%set(gca, 'XScale', 'log')
set(gca, 'FontSize', 16)
xlabel('\it Contrast (Michelson)');
ylabel('\it SVM Accuracy');
title('Psychometric Function Approximations')

sensitivity = 1./thresholdContrasts;

% Visualize relationship between spatial frequency and sensitivity (CSF)
figure(2)
frequencyRange
load compute/ts_CSF_M.mat
frequencyRange
plot( frequencyRange , sensitivity/toDivide,'b.-','MarkerSize',20)
hold on
plot(ts1(:,1),ts1(:,2),'k.-')
plot(ts2(:,1),ts2(:,2),'kx-')
plot(ts3(:,1),ts3(:,2),'ko-')
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gca,'ylim',[0,max(sensitivity)])

if max(frequencyRange) > 2
    set(gca,'xlim',[.1,max(frequencyRange)])
    contrastTicks = [0.1 0.2 0.5 1.0 2.0, max(frequencyRange) ];
    contrastTickLabels = {'.1', '.2', '.5', '1', '2', sprintf('%.0f',max(frequencyRange))};
else
    set(gca,'xlim',[.1,2])
    contrastTicks = [0.1 0.2 0.5 1.0 2.0 ];
    contrastTickLabels = {'.1', '.2', '.5', '1', '2'};
end
    set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 16)
xlabel('\it Spatial Frequency');
ylabel('\it Sensitivity');
title(sprintf('CSF, n = %.2f',nTrialsNum))


%%

if psychFn
    %
    % %at this point, we have a psychometric function of accuracy as a function
    % %of contrast that should sample the threshold well
    
    % if I want to, I can compute an actual psychometric function here:
    
    %pick the first one I chose
    %%
    j = 1;
    numSearchPoints = 5;
    
    nTrials = nTrialsNum/10;
        
    while 1
        a = length(contrastsTotal{1,j});
        if a>10
            break;
        end
        j = j+1;
    end
    
    contrastsTemp = contrastsTotal{1,j};
    
    thresholdContrasts = contrastsTemp(length(contrastsTemp) - ...
        numSearchPoints:length(contrastsTemp));
    
    minC = min(thresholdContrasts);
    maxC = max(thresholdContrasts);
    rangeC = 0.005;%max(thresholdContrasts) - min(thresholdContrasts);
    
    minSample = minC - 2 * rangeC;
    maxSample = maxC + 2 * rangeC;
    
    contrastsToPlot = minSample : rangeC/2 : maxSample
    %%
    if compute == false

        % Size of testing set, based on 10-fold cross-evaluation
        
        % Create presentation display
        presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
        
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
    
    end
    
    contrastsPlot = [];
    accuraciesPlot = [];
    for i = 1:length(contrastsToPlot)
        
        stimParams.contrast = contrastsToPlot(i);
        
        testScene = generateGaborScene(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay);
        
        [acc,time] = getSVMAcc(theMosaic, testScene, nullScene, nTrialsNum);
        
        contrastsPlot = [contrastsPlot, contrastsToPlot(i)];
        accuraciesPlot = [accuraciesPlot, acc];
        
    end
    
    [contrastThreshold,hiResContrasts,hiResPerformance] = getPsychometricFit(contrastsPlot,accuraciesPlot,nTrials);
    %%
    plot(hiResContrasts, hiResPerformance, 'r-', 'LineWidth', 1.5);
    %%
    hold on
    
    plot(contrastsPlot,accuraciesPlot, 'ko', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.7 0.5 0.5], 'MarkerEdgeColor', [0.5 0 0]);
    %%
    set(gca, 'YLim', [0.4 1.05], 'XScale', 'log')
    hold off
    
end