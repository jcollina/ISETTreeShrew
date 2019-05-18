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

angleShift = 20; %degrees, in order to get overlap of 120 degrees

leftConeDensityHuman = zeros(1,singleEyeFOV);

ang = (-singleEyeFOV/2):(singleEyeFOV/2);

for i = 1:length(ang)
        leftConeDensityHuman(i) = (54/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(ang(i)), 'eccentricityUnits', 'deg', 'whichEye', 'left');
end

leftConeDensityHuman(isnan(leftConeDensityHuman)) = min(leftConeDensityHuman);

rightConeDensityHuman = zeros(1,singleEyeFOV);

for i = 1:length(ang)
        rightConeDensityHuman(i) = (54/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(ang(i)), 'eccentricityUnits', 'deg', 'whichEye', 'right');
end

rightConeDensityHuman(isnan(rightConeDensityHuman)) = min(rightConeDensityHuman);

% Now, repeat the process for the left eye
vis_degrees = atand(rel_L./ts_FL) - angleShift + returnShift;

%% Now, to add the two cone density functions together:
leftAngHuman = ang - angleShift;
rightAngHuman = ang + angleShift;

visAngleList = [(leftAngHuman)' ; (rightAngHuman)'];

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

plot(uniqueVisualAnglesHuman,log(totalConeDensityHuman),'r-')

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


