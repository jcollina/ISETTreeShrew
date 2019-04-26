[leftAngTS,rightAngTS,uniqueVisualAnglesTS,leftConeDensityTS,rightConeDensityTS, ...
    tsCD] = getTreeshrewCD('B',false,30);

[leftAngHuman,rightAngHuman,uniqueVisualAnglesHuman,leftConeDensityHuman, ...
    rightConeDensityHuman, humanCD] = getHumanCD('ignoreBinocularOverlap',false);

% Plotting

figure
plot(uniqueVisualAnglesHuman, humanCD,'r-', 'Linewidth', 2)
set(gca, 'YScale', 'log')
hold on
plot(uniqueVisualAnglesTS, tsCD, 'b-', 'Linewidth', 2)

plot(rightAngHuman,rightConeDensityHuman, 'r--')
plot(leftAngHuman,leftConeDensityHuman,'r--')
plot(leftAngTS,leftConeDensityTS,'b--')
plot(rightAngTS,rightConeDensityTS,'b--')
legend('Total Human Cone Density','Total Treeshrew Cone Density')


