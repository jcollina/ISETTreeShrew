
ts_AL = 8; %mm

ts_FL = 5; %mm

fineness = 1; %0/1/2: step size of 1.0/0.1/0.01... visual degrees

angleShift = 58;

returnShift = 0;

%These vectors directly adapted from Muller et al 1989
%ts_fixed_loc = [0, 1.5, 5, 7, 9, 10, 10.5]; %mm
%ts_CD_coarse = [16, 20, 24, 28, 32, 28, 24]; %1000 cones/mm^2

ts_fixed_loc = [0, 1, 3, 6, 9, 10, 12, 14]; %mm
ts_CD_coarse = [16, 20, 24, 28, 32, 32, 28, 24]; %1000 cones/mm^2

% looking at the figures: where is the center of the pupil? 
%   
central_Loc = (max(ts_fixed_loc)-min(ts_fixed_loc))/2; 


rel_R = ts_fixed_loc - central_Loc;
%

maxR = max(abs(rel_R));

a = (maxR/(ts_AL/2));
proj_ecc = (ts_AL/2) * sin(a);
dist = (ts_AL/2) * cos(a);

theor_FL = proj_ecc/tand(88) - dist + (ts_AL/2);

a = (rel_R ./ (ts_AL/2));
proj_ecc = (ts_AL/2) .* sin(a);
dist = (ts_AL/2) .* cos(a);
focal_dist = dist - (ts_AL/2) + theor_FL;
vis_degrees = atand(proj_ecc ./ focal_dist) + angleShift - returnShift;

ts_VD_R = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
ts_CD_R = interp1(vis_degrees,ts_CD_coarse,ts_VD_R,'spline');

%range(vis_degrees)

%%

rel_L = central_Loc - ts_fixed_loc;
%vis_degrees = atand(rel_L/ts_FL);

a = (rel_L ./ (ts_AL/2));
proj_ecc = (ts_AL/2) .* sin(a);
dist = (ts_AL/2) .* cos(a);
focal_dist = dist - (ts_AL/2) + theor_FL;
vis_degrees = atand(proj_ecc ./ focal_dist) - angleShift + returnShift;

ts_VD_L = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
ts_CD_L = interp1(vis_degrees,ts_CD_coarse,ts_VD_L,'spline');

%%
m_L = [ts_VD_L',ts_CD_L'];
m_R = [ts_VD_R',ts_CD_R'];

A = [m_L ; m_R];
A = round(A,fineness);
uniqueCol1 = unique(A(:,1));
%%
% Find out how many times each number occurs.
counts = hist(A(:,1), uniqueCol1);

% Find out how big our output matrix needs to be, and create it.
rows = length(uniqueCol1);
columns = 1 + max(counts);
outputA = zeros(rows, columns);%, 'int32');

% Assign first column
outputA(:,1) = uniqueCol1;
cols234 = A(:,2:end);

% Fill up output rows.
for row = 1 : length(uniqueCol1)
  % Get a unique number from column 1 of A
  thisNumber = uniqueCol1(row);
  rowsWithThisNumber = A(:,1) == thisNumber; % Logical vector.
  theseRows = cols234(rowsWithThisNumber,:)'; % Row vector
  % Assign this row vector to the row from col 2 onwards.
  outputA(row,2:numel(theseRows)+1) = theseRows(:);
end

visualDegrees = outputA(:,1);

totalConeDensity = sum(outputA(:,2:3),2);

figure

plot(visualDegrees,totalConeDensity)

hold on

plot(ts_VD_R,ts_CD_R)

plot(ts_VD_L,ts_CD_L)
