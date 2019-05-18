% Calculations for plotting human and treeshrew cone densities, simplified
% to a single eye

% Which retina (from Muller, 1989) are you interested in analyzing

retinaChoose = 'B';
returnShift = 0;

ts_FL = 5.81;

humanIgnoreOverlap = true;

[rightAngTS,rightConeDensityTS] = getRightTreeshrewCD(retinaChoose);

[rightAngHuman,rightConeDensityHuman] = getRightHumanCD();

% Plotting
figure(1)

plot(rightAngHuman,rightConeDensityHuman, 'r--')
hold on
plot(rightAngTS,rightConeDensityTS,'b--')

set(gca, 'YScale', 'log')
legend('Human Cone Density','Treeshrew Cone Density')
title(['treeshrew ',retinaChoose])
%xlim([0,20])
xlabel(['Visual Angle (',char(176),')'])
ylabel('Combined Cone Density (cones/degree)')