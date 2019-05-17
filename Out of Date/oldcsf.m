function t_BinarySearchCSF(varargin)
%
%% t_BinarySearchCSF (Human/Treeshrew)
% Use a binary SVM and binary search algorithm to create a contrast
% sensitivity function.
%
% Description:
%   Simulating Casagrande's 1984 CSF method (used with
%   treeshrews) in ISETBio. However, the steps can be used for simulating
%   the task on humans optics as well. The key steps are as follows:
%
%   1) Choose the specific spatial frequencies for which you will
% determine the sensitivity to contrast. Currently, ISETBio can only model
% the high-frequency CSF dropoff, not the low-frequency dropoff. Casagrande
% found that the high-frequency dropoff for tree shrews was between 0.75
% and 2 cycles/degree.
%
%   2) For each spatial frequency, determine the contrast level (AKA
% the contrast threshold) that causes the tree shrews to accurately discern
% the contrast 75% of the time.
%
%     * In ISETBio, this step is accomplished by using a support vector
%       machine (SVM) approach to determine the fraction of the time that
%       the algorithm can detect a difference in the mosaic excitation
%       between a null scene and a scene with a contrast. A binary search
%       algorithm is used to identify the specific contrast that leads to
%       75% of the SVM predictions being accurate.
%
%       The binary search works as follows: The spatial frequency of the
%       gabor stimulus is set at the beginning of the search, and is held
%       constant throughout. A minimum and maximum contrast is chosen, with
%       the assumption that the 75% contrast is within the chosen range.
%       The steps involved in the binary search process are well-described
%       in https://en.wikipedia.org/wiki/Binary_search_algorithm.
%
%       Because there is some variability in the SVM results, it is
%       possible for the binary search algorithm to focus on the wrong
%       contrast, without reaching the desired accuracy. For 1000 trials
%       and 10 folds, the SVM results usually have a standard error between
%       1% and 2%. We approach this problem by setting a maximum number of
%       iterations. If that number of iterations is reached without the
%       desired contrast being found, then the code checks for the contrast
%       that led to an accuracy closest to 75%. If 75% falls within a 95%
%       CI for that SVM, then we use that contrast as our "threshold
%       contrast" and move on to the next spatial frequency.
%
%   3) Calculate the contrast sensitivity at each spatial frequency by
% taking the reciprocal of the threshold contrast.
%
%   4) Plot the sensitivity as a function of the spatial frequency.
%
% Optional key/value pairs:
%   csfData - Character array or structure input referencing previously
%       gathered data, if you are specifically interested in plotting that.
%   dataName - If you have a specific choice for the name of the generated
%       data, you can input dataName as a character array.
%   overwrite - If there is already a dataset with the same name, do you
%   want to overwrite it? true/false.

% See also:
%   coneMosaicTreeShrewCreate
%   ls_inferenceTreeShrewBinarySVM
%   t_psychometricFromBinarySearch

% Tools used: getSVMAcc, plotCSF

% History:
%   03/14/19 jsc  Wrote initial version.

%% Initialize workspace and close old figures
%clear; close all;
%ieInit;

%% What do you want to do?
%
% If you're not interested in computing a new CSF, but instead want to
% simply plot prevously gathered data, then input either the name of the
% CSF dataset or the data structure itself into the function. Otherwise,
% leave the input field empty.

%% Parameters

% Do you want to plot Casagrande's data as well?
expData = 'TRUE';

% What species are you interested in simulating? Choose 'treeshrew' or
% 'human'.
species = 'treeshrew';

% How much data do you want to use to train the SVM?
nTrialsNum = 250;

% If you're simulating tree shrews, what do you want the cone density
% of your mosaic to be? This is controlled by changing the minimum cone
% separation.
%   - min cone separation of 6 um     ->      ~ 32,000 cones/mm^2
%   - min cone separation of 7.5 um   ->      ~ 22,000 cones/mm^2
%   - min cone separation of 8.5 um   ->      ~ 16,000 cones/mm^2
cone_spacing = 7.5; %um

% What do you want the standard deviation of the point-spread function
% to be? For humans, this value is approximately 7 um. For tree shrews,
% our preliminary analyses have indicated the value is around 12 um.
psfSigma = 12; %um

% How large do you want the stimulus, and therefore the activated cone
% mosaic, to be? The Casagrande paper used a stimulus of ~14 x 14
% degrees. However, this would take a very long time to run the SVM, so
% we recommend using a size of between 4x4 and 7x7 degrees
sizeDegs = 5; % degrees per side

% What discrete spatial frequencies do you want to find the
% sensitivities for? The Casagrande paper showed a treeshrew CSF
% dropoff in the range of 0.75 - 2 cycles/degree, while for humans it
% is in the range of 1 - 10 cycles/degree.
frequencyRange = 0.75:0.25:2; %cycles per degree

% What range of contrasts do you want to explore for the spatial
% frequencies? Make it wide enough to include accuracies other than
% 100% and 50%!
contrastRange = [.001,.03];

% For the binary search, what range do you want it to reach before
% stopping?
desiredAccRange = [74.7,75.3];

% For the binary search, how many iterations do you want before
% stopping? As appealing as it may be to set a high number here, it is
% important to remember that a mere 10 iterations will lead to 1/1024
% of the original search range- likely past the precision of the SVM.
maxCycles = 10;

% Parse optional input from function
p = inputParser;
p.StructExpand = false;
p.addParameter('csfData', false)
p.addParameter('dataName', char.empty, @ischar)
p.addParameter('overwrite', false, @islogical)
p.parse(varargin{:});
csfData = p.Results.csfData;
dataName = p.Results.dataName;
overwrite = p.Results.overwrite;

% Load a previously generated set of data, if indicated.
switch class(csfData)
    case 'logical'
        compute = true;
    case 'struct'
        compute = false;
    case 'char'
        compute = false;
        csfData = load(csfData);
    otherwise
        error('csfData should either be a data structure of the name of a data structure.')
end

if isstruct(csfData)&&(~isempty(dataName))
    warning('You are loading data instead of generating new data, so no new data will be saved under your chosen name.')
end

% When plotting, it is important to remember that the sensitivity of the
% ideal observer will be much higher than the behaviorally-observed
% sensitivity. This has been shown to be at least a factor of 20. The
% parameter toDivide simply shifts the sensitivity function downward before
% plotting.

toDivide = 20;

if compute
    
    %% Initialize Variables
    %
    % Initialize cell array to store psychometric function data
    contrastsTotal = cell(1,length(frequencyRange));
    accuraciesTotal = cell(1,length(frequencyRange));
    
    % Initialize vector to store threshold contrasts for each spatial
    % frequency
    thresholdContrasts = zeros(1,length(frequencyRange));
    finalAccuracy = zeros(1,length(frequencyRange));
    finalSE = zeros(1,length(frequencyRange));
    
    %% Create a cone mosaic
    %
    % Create a cone mosaic of the appropriate size and cone density. If
    % you're interested in modeling the optics of another species, this is
    % the first place you would make changes.
    switch species
        case 'treeshrew'
            theMosaic = coneMosaicTreeShrewCreate(75, ...
                'fovDegs', sizeDegs, ...
                'customLambda', cone_spacing);
            theOI = oiTreeShrewCreate(...
                'inFocusPSFsigmaMicrons', psfSigma ...
                );
        case 'human'
            theMosaic = coneMosaicHex(7, ...
                'eccBasedConeDensity', true, ...
                'maxGridAdjustmentIterations',200, ...
                'fovDegs', sizeDegs);
            theOI = oiCreate(...
                'inFocusPSFsigmaMicrons', psfSigma ...
                );
        otherwise
            error('species should be treeshrew or human')
    end
    
    %% Main loop
    %
    % Begin looping over the discrete spatial frequencies chosen
    fprintf('Progress Across Frequencies:\n');
    fprintf(['\n' repmat('.',1,length(frequencyRange)) '\n\n']);
    T = tic;
    
    % Matlab generates warnings when you create variables during parallel
    % for loops, so we suppress this specific warning during the loop.
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
    parfor i = 1:length(frequencyRange)
        
        % Create presentation display
        presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
        % Create parameter structure for a low spatial frequency Gabor
        % stimulus
        stimParams = struct(...
            'spatialFrequencyCyclesPerDeg', frequencyRange(i), ...  % changing cycles/deg
            'orientationDegs', 0, ...               % 0 degrees
            'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
            'sizeDegs', sizeDegs, ...               % size
            'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
            'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
            'pixelsAlongWidthDim', [], ...          % pixels- width dimension
            'pixelsAlongHeightDim', [] ...          % pixel- height dimension
            );
        
        % Generate null scene (contrast = 0)
        stimParams.contrast = 0;
        
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
            
            % Generate test stimulus based on current contrast of interest
            currentContrast = mean([maxContrast,minContrast]);
            
            stimParams.contrast = currentContrast;
            
            testScene = generateGaborScene(...
                'stimParams', stimParams,...
                'presentationDisplay', presentationDisplay, ...
                'minimumPixelsPerHalfPeriod', 5);
            
            % Use a binary SVM to get a measure of accuracy in terms of
            % determining which stimulus had gratings
            svm_results = getSVMAcc(theMosaic, theOI, testScene, nullScene, 'nTrialsNum', nTrialsNum);
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
            elseif acc > max(desiredAccRange)
                maxContrast = currentContrast;
            else
                % if reached desired range
                thresholdContrast = currentContrast;
                finalAcctemp = acc;
                finalSEtemp = acc_SE;
                %disp('next :)')
                break;
            end
            
            % If you don't make it, find the closest contrast- if 75%
            % is within the SE of that SVM, take that contrast
            if cycle > maxCycles
                temp = [abs(75-rmmissing(accuracies)) ; rmmissing(accuracies) ; rmmissing(contrasts) ; rmmissing(standardError)]';
                temp = sortrows(temp,1);
                
                tempAcc = temp(1,2);
                thresholdContrast = temp(1,3);
                tempSE = temp(1,4);
                
                if (tempAcc-2*tempSE) < 75 && (tempAcc+2*tempSE) > 75
                    finalAcctemp = tempAcc;
                    finalSEtemp = tempSE;
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
        
        rm(thresholdContrast,finalAcctemp,finalSEtemp)
        
        fprintf('\b|\n');
    end
   
    time = toc(T)/60;
    disp(['The binary search took ', num2str(time), ' minutes.'])
    
    %% Save Data
    
    if isempty(dataName)
        csfDataToSave = [species,sprintf('_csf_%.0f_trials.mat',nTrialsNum)];
    elseif endsWith(dataName,'.mat')
        csfDataToSave = dataName;
    else
        csfDataToSave = strcat(dataName,'.mat');
    end
    
    % Will the data be overwritten?
    if ~overwrite        
        k = 1;
        csfDataToSaveTemp = csfDataToSave;
        while exist(csfDataToSaveTemp, 'file')
            csfDataToSaveTemp = strcat('v',num2str(k),'_',csfDataToSave);
            k = k + 1;
        end
        csfDataToSave = csfDataToSaveTemp;
    end
    
    %create binary search structure
    binaryResults.frequencyRange = frequencyRange;
    binaryResults.contrastsTotal = contrastsTotal;
    binaryResults.accuraciesTotal = accuraciesTotal;
    binaryResults.thresholdContrasts = thresholdContrasts;
    binaryResults.finalAccuracy = finalAccuracy;
    binaryResults.finalSE = finalSE;
    binaryResults.time = time;
    
    save(csfDataToSave,'species','nTrialsNum','psfSigma','theMosaic','theOI','binaryResults')

    csfData = load(csfDataToSave);
end
% Plot csfData
%
% Visualize relationship between contrast and accuracy.

% First, plot the iterations of each binary search in order to make sure
% the SVM results are monotonic. The threshold contrasts are marked by
% asterisks.
plotBinarySearch(csfData.binaryResults)

% If the SVM results are monotonic, it means that the ~75% contrast found
% in the binary search is the only contrast with that accuracy. Now, we can
% plot each threshold contrast as a function of the spatial frequency.

plotCSF(csfData,'toDivide',toDivide,'expData',expData)

%% Functions
    function plotBinarySearch(data)
        
        maxCont = max(cellfun(@(x) max(x),data.contrastsTotal));
        minCont = min(cellfun(@(x) min(x),data.contrastsTotal));
        
        figure()
        hold on
        for j = 1:length(data.frequencyRange)
            tempMat = [data.contrastsTotal{1,j};data.accuraciesTotal{1,j}];
            mat = sortrows(tempMat');
            color = rand(1,3);
            plot(mat(:,1),mat(:,2),'LineWidth',2,'color',color)
            plot(data.thresholdContrasts(j),data.finalAccuracy(j),'*','MarkerSize',10,'color',color)
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

    function plotCSF(data,varargin)
        
        p = inputParser;
        p.addParameter('toDivide', 1 , @isnumeric);
        p.addParameter('expData','FALSE', @ischar);
        
        % Parse input
        p.parse(varargin{:});
        
        toDivide = p.Results.toDivide;
        expData = p.Results.expData;
        
        % Shifting the sensitivity downwards by a factor of 'toDivide', if desired
        sensitivity = (1./data.binaryResults.thresholdContrasts)./toDivide;
        
        plotTitle = sprintf('CSF \n %.0f trials, sigma_{PSF} = %.0f, mosaic size = %.0f',data.nTrialsNum,data.psfSigma,round(data.theMosaic.fov));
        
        figure()
        
        plot(data.binaryResults.frequencyRange , sensitivity,'b.-','MarkerSize',20)
        
        % If you want to plot data from Casagrande 1984:
        if expData
            hold on
            load ts_CSF_M.mat ts_CSF_M
            
            ts1 = ts_CSF_M{1};
            ts2 = ts_CSF_M{2};
            ts3 = ts_CSF_M{3};
            
            plot(ts1(:,1),ts1(:,2),'k.-')
            plot(ts2(:,1),ts2(:,2),'kx-')
            plot(ts3(:,1),ts3(:,2),'ko-')
            
            legend('ISETTreeShrew CSF','Casagrande TS0','Casagrande TS1','Casagrande TS3','Location','northwest')
            set(gca,'ylim',[0,max([sensitivity,ts3(:,2)'])])
        else    
            set(gca,'ylim',[0,max(sensitivity)])
        end
        
        if max(data.binaryResults.frequencyRange) > 2
            set(gca,'xlim',[.1,max(data.binaryResults.frequencyRange)])
            contrastTicks = [0.1 0.2 0.5 1.0 2.0, max(data.binaryResults.frequencyRange) ];
            contrastTickLabels = {'.1', '.2', '.5', '1', '2', sprintf('%.0f',max(data.binaryResults.frequencyRange))};
        else
            set(gca,'xlim',[.1,2])
            contrastTicks = [0.1 0.2 0.5 1.0 2.0 ];
            contrastTickLabels = {'.1', '.2', '.5', '1', '2'};
        end
        
        set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        set(gca, 'FontSize', 16)
        xlabel('\it Spatial Frequency (Cycles/Degree)');
        ylabel('\it Contrast Sensitivity');
        title(plotTitle)
    end
end