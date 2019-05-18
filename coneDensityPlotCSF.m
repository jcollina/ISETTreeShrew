function coneDensityPlotCSF

toDivide = 20;


datareg = load('csf_1000_trials_psf_12_size_5.mat');
datamin = load('2cd_min_csf_1000_trials_psf_12_size_5.mat');
datamax = load('2cd_max_csf_1000_trials_psf_12_size_5.mat');

sensitivityreg = 1./datareg.thresholdContrasts;
sensitivitymin = 1./datamin.thresholdContrasts;
sensitivitymax = 1./datamax.thresholdContrasts;

% Visualize relationship between spatial frequency and sensitivity (CSF)

%load ts_CSF_M.mat ts_CSF_M

ts1 = ts_CSF_M{1};
ts2 = ts_CSF_M{2};
ts3 = ts_CSF_M{3};

%plotTitle = sprintf('CSF \n %.0f trials, sigma_{PSF} = %.0f, mosaic size = %.0f',data.nTrialsNum,data.psfSigma,data.sizeDegs);

sens = [sensitivitymax, sensitivitymin, sensitivityreg]./toDivide;

figure()
hold on
set(gca,'ylim',[0,2*max(sens)],'YScale','log')
set(gca,'xlim',[.75,2],'XScale', 'log')
contrastTicks = [0.1 0.2 0.5 1.0 2.0 ];
contrastTickLabels = {'.1', '.2', '.5', '1', '2'};

set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels);
plot( data.frequencyRange , sensitivitymin/toDivide,'r.-','MarkerSize',20)
plot( data.frequencyRange , sensitivityreg/toDivide,'b.-','MarkerSize',20)
plot( data.frequencyRange , sensitivitymax/toDivide,'g.-','MarkerSize',20)
plot(ts1(:,1),ts1(:,2),'k.-')
plot(ts2(:,1),ts2(:,2),'kx-')
plot(ts3(:,1),ts3(:,2),'ko-')

legend('12,000 cones/mm^2','22,000 cones/mm^2','22,000 cones/mm^2','TS0','TS1','TS3')


xlabel('\it Spatial Frequency');
ylabel('\it Sensitivity');
%title(plotTitle)

end