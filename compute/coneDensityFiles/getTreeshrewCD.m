
function [leftAngTS,rightAngTS,uniqueVisualAnglesTS,leftConeDensityTS,rightConeDensityTS, tsCD] = getTreeshrewCD(ts_choose,useTheorFL,returnShift)

% In Muller's 1989 TS Topography Paper, Figure 5 shows the cone density
% maps for 3 separate treeshrew retinas. Out of retinas 'A', 'B' and 'C',
% which do you want to plot the cone density as a function of visual angle?
%ts_choose = 'A';

% What do you want to set as the treeshrew axial length? Norton et al, 1991
% found an AL of approximately 7.8 mm.


% What do you want to set as the treeshrew axial length? Norton et al, 1991
% found an AL of approximately 7.8 mm.
ts_AL = 7.8; %mm 

% We make the simplifying assumption that the focal length and the
% posterior nodal distance (PND, mm from posterior part of lens to the
% retina) are the same.

%do you want to use data that has been gathered about the treeshrew focal
% length, and see what field of view is returned? Or do you want to find
% the focal length that explains a hypothesis about the field of view?
useTheorFL = true;

% If you just want to use the PND used in Jarvis et al, 2003, and see what
% FOV emerges:

ts_FL_exp = 5.81; %mm 

%Question: do we want to use the anterior or posterior focal lengths?
% so: 4.35 mm or 5.81 mm?

% If you want to do this, what is your the field of binocular overlap you
% want to achieve? This value is for the total visual field, not just of
% one eye.
binocularField = 60; %degrees

% 0/1/2: step size of 1.0/0.1/0.01... 1 -> step of 0.1 is usally enough.
fineness = 1; 

% How do you want to interpolate the experimental data? Muller has gradient
% lines, but we want to fill in the rest. I recommend 'pchip' if you don't
% want to extrapolate past the experimental data, 'spline' if you want a
% smoother cone density curve over eccentricity.
interpolationMethod = 'spline';

% What is the angular displacement for the treeshrew eyes? Barnes Jannuzi
% measured the angle of the skull as 58 degrees.
angleShift = 58; %degrees

% If you want to, what is the amount of forward rotation they have? If you
% think they can't rotate their eyes, this is zero.
%returnShift = 20; %degrees

% These vectors are directly adapted from Muller et al 1989. Locations in
% mm are with respect to left side of the slice.
switch ts_choose
    case 'A'
        exp_eye = 'right';
        ts_fixed_loc = [0, 1, 3, 6, 9, 10, 12, 14, 15]; %mm 
        ts_CD_coarse = (range(ts_fixed_loc)/2*1000/(angleShift+binocularField/2)).*[16, 20, 24, 28, 32, 32, 28, 24, 20]; %cones/degree
    case 'B'
        exp_eye = 'left';
        ts_fixed_loc = [0, 0.5, 2.5, 3, 3.5, 5,7, 8, 10, 11.5, 13, 14, 14.25, 14.5, 15]; %mm 
        ts_CD_coarse = (range(ts_fixed_loc)/2*1000/(angleShift+binocularField/2)).*[16, 20, 24, 28, 32, 32, 32, 28, 28, 32, 32, 28, 24, 20, 16]; %cones/degree
    case 'C'
        exp_eye = 'left';
        ts_fixed_loc = [0,0.5,1,2,3,4.5,8,9,12,13,14,15]; %mm 
        ts_CD_coarse = (range(ts_fixed_loc)/2*1000/(angleShift+binocularField/2)).*[16, 20, 24, 28, 32, 36, 36, 32, 28, 24, 20, 16]; %cones/degree
end

% Looking at the figures: where is the center of the pupil? Not explicitly
% stated in paper, but seems to be in center of slice.
central_Loc = (max(ts_fixed_loc)-min(ts_fixed_loc))/2;

% Now, we want the locations to be with respect to the central location. We
% want the left side of the visual field to be negative, and the right side
% to be positive. We therefore need to discriminate between the left and
% right eye.

switch exp_eye
    case 'right'
        rel_R = ts_fixed_loc - central_Loc;
        rel_L = central_Loc - ts_fixed_loc;
    case 'left'
        rel_L = ts_fixed_loc - central_Loc;
        rel_R = central_Loc - ts_fixed_loc;
end

% If desired, we use the basic parameters to calculate the focal length
% that would result in the desired binocular field and total FOV.

if useTheorFL == true
    
    maxEcc = (angleShift-returnShift) + binocularField/2;
    maxR = max(abs(rel_R));
    
    a = (maxR/(ts_AL/2));
    proj_ecc = (ts_AL/2) * sin(a);
    dist = (ts_AL/2) * cos(a);
    ts_FL = proj_ecc/tand(maxEcc) - dist + (ts_AL/2);
    
else
    ts_FL = ts_FL_exp;
end

a = (rel_R ./ (ts_AL/2));
proj_ecc = (ts_AL/2) .* sin(a);
dist = (ts_AL/2) .* cos(a);
focal_dist = dist - (ts_AL/2) + ts_FL;
vis_degrees = atand(proj_ecc ./ focal_dist) + angleShift - returnShift;

% Now, use the interp1 function to create a discrete vector of cone density
% as a function of visual angle.

% Visual degrees:
rightAngTS = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
% Cone density:
rightConeDensityTS = interp1(vis_degrees,ts_CD_coarse,rightAngTS, interpolationMethod);

% Now, repeat the process for the left eye
vis_degrees = atand(-1*proj_ecc ./ focal_dist) - angleShift + returnShift;

leftAngTS = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
leftConeDensityTS = interp1(vis_degrees,ts_CD_coarse,leftAngTS, interpolationMethod);

% % Now, to add the two cone density functions together:

visAngleList = round([leftAngTS' ; rightAngTS'],fineness);

coneDensityList = [leftConeDensityTS' ; rightConeDensityTS'];

uniqueVisualAnglesTS = unique(visAngleList);

%% Combining information from both eyes
% Find out how many times each number occurs. Code taken from:
%   https://www.mathworks.com/matlabcentral/answers/152536-how-to-find-the-
%   ...rows-with-the-same-values-and-to-merge-those#comment_234046

% Find out how big our output matrix needs to be, and create it.
combinedMatrix = zeros(length(uniqueVisualAnglesTS), 3);

% We want the first column to be the unique visual angles
combinedMatrix(:,1) = uniqueVisualAnglesTS;

% Loop through visual angles and sum associated cone densities
for i = 1 : length(uniqueVisualAnglesTS)
    angleTemp = uniqueVisualAnglesTS(i);
    tempIndex = visAngleList == angleTemp;
    coneDensitiesTemp = coneDensityList(tempIndex,:)'; % Row vector
    combinedMatrix(i,2:numel(coneDensitiesTemp)+1) = coneDensitiesTemp(:);
end

tsCD = sum(combinedMatrix(:,2:3),2);
end