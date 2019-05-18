%% t_humanEyeGeometry Illustrate effect of various geometric properties of
% the human eye on cone density as a function of external visual angle.
% 

% 
% See also: t_eyeGeometry2D
% 

% History:
%   04/19/19 jsc  Wrote initial version. 

%%

plotTS = 0;

singleEyeFOV = 160;

binocularField = 120;

peripheralField = 40;

AL = 24;
retMM = (0.72*pi*AL);

leftConeDensityHuman = zeros(1,singleEyeFOV);

pos = (-retMM/2):0.1:(retMM/2);

for i = 1:length(pos)
        leftConeDensityHuman(i) = (retMM/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(pos(i)), 'eccentricityUnits', 'mm', 'whichEye', 'left');
end

leftConeDensityHuman(isnan(leftConeDensityHuman)) = min(leftConeDensityHuman);
%%
    maxEcc = (peripheralField+(binocularField/2))/2;
    maxMM = max(abs(pos));
    %%
    a = (maxMM/(AL/2));
    proj_ecc = (AL/2) * sin(a);
    dist = (AL/2) * cos(a);
    FL = proj_ecc/tand(maxEcc) - dist + (AL/2);
%%
a = (pos ./ (AL/2));
proj_ecc = (AL/2) .* sin(a);
dist = (AL/2) .* cos(a);
focal_dist = dist - (AL/2) + FL;
%%
leftAngleHumanTemp = atand(-1.*proj_ecc ./ focal_dist) + angleShift - returnShift;

leftAngleHuman = min(round(leftAngleHumanTemp)):(10^(-1)):max(round(leftAngleHumanTemp));

leftConeDensityHuman = interp1(leftAngleHumanTemp,leftConeDensityHuman,leftAngleHuman, interpolationMethod);

rightConeDensityHuman = zeros(1,singleEyeFOV);

for i = 1:length(pos)
        rightConeDensityHuman(i) = (retMM/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(pos(i)), 'eccentricityUnits', 'deg', 'whichEye', 'right');
end
rightConeDensityHuman(isnan(rightConeDensityHuman)) = min(rightConeDensityHuman);
%%
rightAngleHumanTemp = atand(proj_ecc ./ focal_dist) - angleShift + returnShift;
rightAngleHuman = min(round(rightAngleHumanTemp)):(10^(-1)):max(round(rightAngleHumanTemp));

rightConeDensityHuman = interp1(rightAngleHumanTemp,rightConeDensityHuman,rightAngleHuman, interpolationMethod);
%%

%% Now, to add the two cone density functions together:
leftAngHuman = leftAngleHuman;
rightAngHuman = rightAngleHuman;

visAngleList = round([(leftAngHuman)' ; (rightAngHuman)'],2);

coneDensityList = [leftConeDensityHuman' ; rightConeDensityHuman'];

uniqueVisualAnglesHuman = unique(visAngleList);

%% Combining information from both eyes
% Find out how many times each number occurs. Code taken from:
%   https://www.mathworks.com/matlabcentral/answers/152536-how-to-find-the-
%   ...rows-with-the-same-values-and-to-merge-those#comment_234046

% Find out how big our output matrix needs to be, and create it.
combinedMatrix = zeros(length(uniqueVisualAnglesHuman), 3);

% We want the first column to be the unique visual angles
combinedMatrix(:,1) = uniqueVisualAnglesHuman;

% Loop through visual angles and sum associated cone densities
for i = 1 : length(uniqueVisualAnglesHuman)
    angleTemp = uniqueVisualAnglesHuman(i);
    tempIndex = visAngleList == angleTemp;
    coneDensitiesTemp = coneDensityList(tempIndex,:)'; % Row vector
    combinedMatrix(i,2:numel(coneDensitiesTemp)+1) = coneDensitiesTemp(:);
end

totalConeDensityHuman = sum(combinedMatrix(:,2:3),2);
dif = combinedMatrix(:,3) - combinedMatrix(:,2);
%% Plotting

figure
%plot(rightAngHuman,log(rightConeDensityHuman),'r--')
plot(uniqueVisualAnglesHuman,log(totalConeDensityHuman),'r-')
%{
hold on
if plotTS
plot(uniqueVisualAnglesTS,log(totalConeDensityTS),'b-')
end
plot(leftAngHuman,log(leftConeDensityHuman),'r--')
plot(rightAngHuman,log(rightConeDensityHuman),'r--')
if plotTS   
plot(leftAngTS,log(leftConeDensityTS),'b--')
plot(rightAngTS,log(rightConeDensityTS),'b--')
legend('Total Human Cone Density','Total Treeshrew Cone Density')
else
    legend('Total Human Cone Density','Individual Eye Cone Density')
end
%}

