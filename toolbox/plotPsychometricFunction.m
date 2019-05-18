function plotPsychometricFunction(featureSamples,performance,xLabel,varargin)
% Plot the psychometric function for a given binary search and weibull fit
%
% Syntax:
%   plotPsychometricFunction(featureSamples,performance,'Contrasts (Michelson)','fitResults','psychFitResults')
%
% Inputs:
%   xLabel
%   featureSamples - numeric array, levels of the feature that were sampled
%   in the search
%   performance - numeric array, SVM results associated with the feature
%   xLabel
  
%
% Optional key/value pairs:
%   fitResults - Matlab structure detailing fit of this data with fields
%   'hiResContrast', 'hiResPerformance' and 'contrastThreshold'
%   svmError - error in SVM accuracy, if included will be displayed as
%   error bars

% See also:

% History:
%   05/14/19 jsc  Wrote initial version.

p = inputParser;
p.addParameter('svmError', 0 ,@isnumeric);
p.addParameter('fitResults',struct,@isstruct);
p.parse(varargin{:});
fitResults = p.Results.fitResults;
svmError = p.Results.svmError;

if ismember('fitResults',p.UsingDefaults)
    fit = false;
else
    fit = true;
end

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
    line([min(featureSamples),contrastThreshold],[75,75])
    plot(featureSamples,performance, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
else
    plot(featureSamples,performance, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
    hold on
end

if ~ismember('svmError',p.UsingDefaults)
    errorbar(featureSamples,performance,thresholdSample.acc_SE)
end

%{
if max(featureSamples) > .031
    set(gca,'xlim',[min(featureSamples),max(featureSamples)])
    contrastTicks = [0.005,0.01 0.02 0.03, max(featureSamples) ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', int2str(round(maxContrastsToPlot))};
else
    set(gca,'xlim',[min(featureSamples),.03])
    contrastTicks = [0.005 0.01 0.02 0.03];
    contrastTickLabels = {'0.005', '.01',  '.02',  '.03'};
end
%}
set(gca,'xlim',[min(featureSamples),max(featureSamples)])
%set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);
set(gca, 'YLim', [40 105], 'XScale', 'log')
set(gca, 'FontSize', 16)
xlabel(xLabel);
ylabel('Performance (% Detected)');
title('Individual Psychometric Function')
hold off
end