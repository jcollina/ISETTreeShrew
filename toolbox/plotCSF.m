function plotCSF(data,varargin)
% Plot the contrast sensitivity function for a given simulated visual
% system front-end.
%
% Syntax:
%   plotCSF(data,'toDivide','20','expData','TRUE')
%
% Inputs:
%   data        Matlab structure with fields 'thresholdContrasts',
%               'frequencyRange', 'nTrialsNum', 'psfSigma' and 'sizeDegs'
%
% Optional key/value pairs:
%   toDivide - Factor to shift sensitivity function downwards, if desired
%   expData - Plot Casagrande's experimental data? 'TRUE' or 'FALSE'

% See also:
%   t_BinarySearchCSF

% History:
%   04/08/19 jsc  Wrote initial version.
%   04/28/19 jsc  Formatted.

p = inputParser;
p.addParameter('toDivide', 1 , @isnumeric);
p.addParameter('expData', 'FALSE' , @ischar);

% Parse input
p.parse(varargin{:});

toDivide = p.Results.toDivide;
expData = p.Results.expData;

% Shifting the sensitivity downwards by a factor of 'toDivide', if desired
sensitivity = (1./data.thresholdContrasts)./toDivide;

plotTitle = sprintf('CSF \n %.0f trials, sigma_{PSF} = %.0f, mosaic size = %.0f',data.nTrialsNum,data.psfSigma,data.sizeDegs);

figure()

plot(data.frequencyRange , sensitivity/toDivide,'b.-','MarkerSize',20)

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
end

set(gca,'ylim',[0,max(sensitivity/toDivide)])

if max(data.frequencyRange) > 2
    set(gca,'xlim',[.1,max(data.frequencyRange)])
    contrastTicks = [0.1 0.2 0.5 1.0 2.0, max(data.frequencyRange) ];
    contrastTickLabels = {'.1', '.2', '.5', '1', '2', sprintf('%.0f',max(data.frequencyRange))};
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
