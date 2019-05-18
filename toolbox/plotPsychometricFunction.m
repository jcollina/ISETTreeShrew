function plotPsychometricFunction(thresholdSample,varargin)
% Plot the psychometric function for a given binary search and weibull fit
%
% Syntax:
%   plotPsychometricFunction(thresholdSample,'fitResults','psychFitResults')
%
% Inputs:
%   thresholdSample - Matlab structure detailing SVM accuracy as a function
%   of stimulus contrast for a set stimulus spatial frequency. With fields
%   'contrasts', 'accuracies', 'acc_SE' and 'frequency'.
  
%
% Optional key/value pairs:
%   fitResults - Matlab structure detailing fit of this data with fields
%   'hiResContrast', 'hiResPerformance' and 'contrastThreshold'

% See also:

% History:
%   05/14/19 jsc  Wrote initial version.

p = inputParser;
p.addParameter('fitResults',struct,@isstruct);
p.parse(varargin{:});
fitResults = p.Results.fitResults;

if ismember('fitResults',p.UsingDefaults)
    fit = false;
else
    fit = true;
end

contrastsToPlot = thresholdSample.contrasts;
accuraciesToPlot = thresholdSample.accuracies;
spatFreq = thresholdSample.frequency;

if fit
    hiResContrasts = fitResults.hiResContrasts;
    hiResPerformance = fitResults.hiResPerformance;
    contrastThreshold = fitResults.contrastThreshold;
end

figure()
if fit
    plot(hiResContrasts, hiResPerformance, 'r-', 'LineWidth', 1.5);
    hold on
    line([contrastThreshold,contrastThreshold],[40,75])
    line([min(contrastsToPlot),contrastThreshold],[75,75])
    plot(contrastsToPlot,accuraciesToPlot, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
else
    plot(contrastsToPlot,accuraciesToPlot, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
    hold on
end

if isfield(thresholdSample,'acc_SE')
    errorbar(contrastsToPlot,accuraciesToPlot,thresholdSample.acc_SE)
end

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
xlabel('Contrast (Michelson)');
ylabel('% SVM Accuracy');
title(sprintf('Individual Psychometric Function for\nSpatial Frequency of %.0f cpd',spatFreq))
hold off
end