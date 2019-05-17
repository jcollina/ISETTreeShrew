function plotBinarySearch(binaryResults)
% Plot the contrast sensitivity function for a given simulated visual
% system front-end.
%
% Syntax:
%   plotBinarySearch(binaryResults)
%
% Inputs:
%   binaryResults - Matlab structure with fields 'frequencyRange',
%   'thresholdContrasts', 'finalAccuracy', 'contrastsTotal',
%   'accuraciesTotal', and 'finalSE'

% See also:
%   t_BinarySearchCSF

% History:
%   04/13/19 jsc  Wrote initial version.

% Determine the total contrast range
maxCont = max(cellfun(@(x) max(x),binaryResults.contrastsTotal));
minCont = min(cellfun(@(x) min(x),binaryResults.contrastsTotal));

figure()
hold on

% Plot the individual search points and their corresponding SVM accuracy
for j = 1:length(binaryResults.frequencyRange)
    tempMat = [binaryResults.contrastsTotal{1,j};binaryResults.accuraciesTotal{1,j}];
    mat = sortrows(tempMat');
    color = rand(1,3);
    plot(mat(:,1),mat(:,2),'LineWidth',2,'color',color)
    
    % Use the SE from the SVM results as an error bar for the final
    % accuracy
    errorbar(binaryResults.thresholdContrasts(j),binaryResults.finalAccuracy(j), binaryResults.finalSE(j),'k*','MarkerSize',10)%,'color',color)
end

% Expand the x-axis if the max contrast is higher than 0.03
if maxCont > .031
    set(gca,'xlim',[minCont,maxCont])
    contrastTicks = [0.005,0.01 0.02 0.03, maxCont ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', sprintf('%.0f',maxCont)};
else
    set(gca,'xlim',[minCont,.03])
    contrastTicks = [0.005 0.01 0.02 0.03];
    contrastTickLabels = {'0.005', '.01',  '.02',  '.03'};
end

set(gca, ...
    'XTick', contrastTicks, ...
    'XTickLabel', contrastTickLabels, ...
    'XScale', 'log', ...
    'FontSize', 16 ...
    );

xlabel('\it Contrast (Michelson)');
ylabel('\it SVM Accuracy');
title(['Psychometric Function Approximations' newline 'Using a Binary Search']);

hold off
end