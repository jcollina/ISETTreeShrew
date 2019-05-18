function [rightAngHuman,rightConeDensityHuman] = getRightHumanCD(~)
% 
% Use the coneDensityReadData function to get cone Density as a function of
% visual angle for the right human eye

% Parameters:
singleEyeFOV = 160;

focalLength = 17;
 
% Calculate Cone Density
ang = (-singleEyeFOV/2):(singleEyeFOV/2);
rightConeDensityHuman = zeros(1,singleEyeFOV);
for i = 1:length(ang)
        rightConeDensityHuman(i) = ((focalLength*3.14)/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(ang(i)), 'eccentricityUnits', 'deg', 'whichEye', 'right');
end

rightAngHuman = ang;
rightConeDensityHuman(isnan(rightConeDensityHuman)) = min(rightConeDensityHuman);
end
