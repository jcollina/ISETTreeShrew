function plotCSF(data,toDivide)

%
sensitivity = 1./data.thresholdContrasts;

% Visualize relationship between spatial frequency and sensitivity (CSF)

load ts_CSF_M.mat ts_CSF_M

ts1 = ts_CSF_M{1};
ts2 = ts_CSF_M{2};
ts3 = ts_CSF_M{3};

plotTitle = sprintf('CSF \n %.0f trials, sigma_{PSF} = %.0f, mosaic size = %.0f',data.nTrialsNum,data.psfSigma,data.sizeDegs);

figure()

plot( data.frequencyRange , sensitivity/toDivide,'b.-','MarkerSize',20)
hold on
plot(ts1(:,1),ts1(:,2),'k.-')
plot(ts2(:,1),ts2(:,2),'kx-')
plot(ts3(:,1),ts3(:,2),'ko-')

set(gca,'ylim',[0,max(sensitivity)])

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
xlabel('\it Spatial Frequency');
ylabel('\it Sensitivity');
title(plotTitle)

end
