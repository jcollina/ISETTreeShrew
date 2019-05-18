retinaChoose = 'B';
returnShift = 0;

%Question: do we want to use the anterior or posterior focal lengths?
% so: 4.35 mm or 5.81 mm?
ts_FL = 5.81;

humanIgnoreOverlap = true;

[leftAngTS,rightAngTS,uniqueVisualAnglesTS,leftConeDensityTS,rightConeDensityTS, ...
    tsCD] = getTreeshrewCD(retinaChoose,ts_FL,returnShift);

[leftAngHuman,rightAngHuman,uniqueVisualAnglesHuman,leftConeDensityHuman, ...
    rightConeDensityHuman, humanCD] = getHumanCD('ignoreBinocularOverlap',humanIgnoreOverlap);

% Plotting
figure(1)
%plot(uniqueVisualAnglesHuman, humanCD,'r-', 'Linewidth', 2)
set(gca, 'YScale', 'log')
hold on
%plot(uniqueVisualAnglesTS, tsCD, 'b-', 'Linewidth', 2)

plot(rightAngHuman,rightConeDensityHuman, 'r--')
hold on
%plot(leftAngHuman,leftConeDensityHuman,'r--')
%plot(leftAngTS,leftConeDensityTS,'b--')
plot(rightAngTS,rightConeDensityTS,'b--')

legend('Total Human Cone Density','Total Treeshrew Cone Density')
title(sprintf('treeshrew %s with eyes shifted %.0f degrees forward', retinaChoose,returnShift))
xlim([0,20])
xlabel(['Visual Angle (',char(176),')'])
ylabel('Combined Cone Density (cones/degree)')