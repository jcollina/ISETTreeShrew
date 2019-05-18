function [leftAngHuman,rightAngHuman,uniqueVisualAnglesHuman,leftConeDensityHuman,rightConeDensityHuman, humanCD] = getHumanCD(varargin)
% Illustrate effect of various geometric properties of
% the human eye on cone density as a function of external visual angle.
% 

p = inputParser;
p.addParameter('ignoreBinocularOverlap',true, @islogical);
p.addParameter('singleEyeFOV', 160, @isnumeric);
% Parse input
p.parse(varargin{:});

ignoreBinocularOverlap = p.Results.ignoreBinocularOverlap;
singleEyeFOV = p.Results.singleEyeFOV;

if ignoreBinocularOverlap
    angleShift = 0;
else
    angleShift = 20; %degrees, in order to get overlap of 120 degrees
end

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
humanCD = sum(combinedMatrix(:,2:3),2);
dif = combinedMatrix(:,3) - combinedMatrix(:,2);

end
