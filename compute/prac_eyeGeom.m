ts_FL = 4.35; %mm
ts_CD_coarse = [16,20,24,28,32,28,24]; %1000 cones/mm^2

ts_R = [0,1.5,5,7,9,10,10.5]; %mm
central_Loc = 1.2;
rel_R = ts_R - central_Loc;
vis_degrees = atand(rel_R/ts_FL);

ts_VD_R = min(vis_degrees):0.1:max(vis_degrees);
ts_CD_R = interp1(vis_degrees,ts_CD_coarse,ts_VD_R,'spline');

ts_L = [0,-1.5,-5,-7,-9,-10,-10.5]; %mm
central_Loc = -1.2;
rel_L = ts_L - central_Loc;
vis_degrees = atand(rel_L/ts_FL);

ts_VD_L = min(vis_degrees):0.1:max(vis_degrees);
ts_CD_L = interp1(vis_degrees,ts_CD_coarse,ts_VD_L,'spline');

%%

figure
plot(ts_VD_L,ts_CD_L,'.')
hold on
plot(ts_VD_R,ts_CD_R,'.')
%%
m_L = [ts_VD_L',ts_CD_L'];
m_R = [ts_VD_R',ts_CD_R'];

[m_L,m_R]