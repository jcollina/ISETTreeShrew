function plotPsychometricFunction(thresholdSample,fitResults)

contrastsToPlot = thresholdSample.contrastsTotal;
accuraciesToPlot = thresholdSample.accuraciesTotal;
acc_SE = thresholdSample.totalSE;
spatFreq = thresholdSample.frequency;

hiResContrasts = fitResults.hiResContrasts;
hiResPerformance = fitResults.hiResPerformance;
contrastThreshold = fitResults.contrastThreshold;


        figure()
        plot(hiResContrasts, hiResPerformance, 'r-', 'LineWidth', 1.5);
        hold on
        errorbar(contrastsToPlot,accuraciesToPlot,acc_SE)
        plot(contrastsToPlot,accuraciesToPlot, 'ko', 'MarkerSize', 8, ...
            'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
        line([contrastThreshold,contrastThreshold],[40,75])
        line([min(contrastsToPlot),contrastThreshold],[75,75])
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