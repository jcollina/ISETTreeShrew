%% t_simpleTSEyeGeometry 
% Illustrate effect of various geometric properties of the treeshrew eye on
% cone density as a function of external visual angle.
% 

% 
% See also: t_eyeGeometry2D
% 

% In Muller's 1989 TS Topography Paper, Figure 5 shows the cone density
% maps for 3 separate treeshrew retinas. Out of retinas 'A', 'B' and 'C',
% which do you want to plot the cone density as a function of visual angle?
ts_choose = 'A';

% What do you want to set as the treeshrew axial length? Norton et al, 1991
% found an AL of approximately 7.8 mm.
ts_AL = 7.8; %mm 

% 0/1/2: step size of 1.0/0.1/0.01... 1 -> step of 0.1 is usally enough.
fineness = 1; 

% How do you want to interpolate the experimental data? Muller has gradient
% lines, but we want to fill in the rest. I recommend 'pchip' if you don't
% want to extrapolate past the experimental data, 'spline' if you want a
% smoother cone density curve over eccentricity.
interpolationMethod = 'pchip';

% What is the angular displacement for the treeshrew eyes? Barnes Jannuzi
% measured the angle of the skull as 58 degrees.
angleShift = 58; %degrees

% If you want to, what is the amount of forward rotation they have? If you
% think they can't rotate their eyes, this is zero.
returnShift = 40;%40; %degrees

% These vectors are directly adapted from Muller et al 1989. Locations in
% mm are with respect to left side of the slice.

switch ts_choose
    case 'A'
        exp_eye = 'right';
        ts_fixed_loc = [0, 1, 3, 6, 9, 10, 12, 14, 15]; %mm 
        ts_CD_coarse = (15*1000/176).*[16, 20, 24, 28, 32, 32, 28, 24, 20]; %cones/degree
    case 'B'
        exp_eye = 'left';
        ts_fixed_loc = [0, 0.5, 2.5, 3, 3.5, 5,7, 8, 10, 11.5, 13, 14, 14.25, 14.5, 15]; %mm 
        ts_CD_coarse = (15*1000/176).*[16, 20, 24, 28, 32, 32, 32, 28, 28, 32, 32, 28, 24, 20, 16]; %cones/degree
    case 'C'
        exp_eye = 'left';
        ts_fixed_loc = [0,0.5,1,2,3,4.5,8,9,12,13,14,15]; %mm 
        ts_CD_coarse = (15*1000/176).*[16, 20, 24, 28, 32, 36, 36, 32, 28, 24, 20, 16]; %cones/degree
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
        rel_R = (ts_fixed_loc - central_Loc);
        rel_L = -1.*rel_R;
    case 'left'
        rel_L = (ts_fixed_loc - central_Loc);
        rel_R = -1.*rel_L;
end

vis_degrees = atand(rel_R./ts_FL) + angleShift - returnShift;

% Now, use the interp1 function to create a discrete vector of cone density
% as a function of visual angle.

% Visual degrees:
ts_VD_R = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
% Cone density:
ts_CD_R = interp1(vis_degrees,ts_CD_coarse,ts_VD_R, interpolationMethod);

% Now, repeat the process for the left eye
vis_degrees = atand(rel_L./ts_FL) - angleShift + returnShift;

ts_VD_L = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
ts_CD_L = interp1(vis_degrees,ts_CD_coarse,ts_VD_L, interpolationMethod);

% % Now, to add the two cone density functions together:

visAngleList = round([ts_VD_L' ; ts_VD_R'],fineness);

coneDensityList = [ts_CD_L' ; ts_CD_R'];

uniqueVisualAngles = unique(visAngleList);

%% Combining information from both eyes
% Find out how many times each number occurs. Code taken from:
%   https://www.mathworks.com/matlabcentral/answers/152536-how-to-find-the-
%   ...rows-with-the-same-values-and-to-merge-those#comment_234046

% Find out how big our output matrix needs to be, and create it.
combinedMatrix = zeros(length(uniqueVisualAngles), 3);

% We want the first column to be the unique visual angles
combinedMatrix(:,1) = uniqueVisualAngles;

% Loop through visual angles and sum associated cone densities
for i = 1 : length(uniqueVisualAngles)
    angleTemp = uniqueVisualAngles(i);
    tempIndex = visAngleList == angleTemp;
    coneDensitiesTemp = coneDensityList(tempIndex,:)'; % Row vector
    combinedMatrix(i,2:numel(coneDensitiesTemp)+1) = coneDensitiesTemp(:);
end

totalConeDensity = sum(combinedMatrix(:,2:3),2);

%% Plotting

figure

plot(uniqueVisualAngles,totalConeDensity)
hold on
plot(ts_VD_L,ts_CD_L)
plot(ts_VD_R,ts_CD_R)

if returnShift == 0
    title(sprintf('Cone Density for treeshrew % s',ts_choose))
else
    title(sprintf('Cone Density for Treeshrew % s \n with Eye Rotation of % .0f% s',ts_choose,returnShift,char(176)))
end

legend('Total Cone Density','Left Eye Cone Density','Right Eye Cone Density')
xlabel(['Visual Angle (',char(176),')'])
ylabel('Combined Cone Density (cones/degree)')