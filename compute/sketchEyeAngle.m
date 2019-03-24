% questions to ask:
% 1) what's known about the specific drop-off in cone density? It
% definitely doesn't seem spherically symmetric, but is there at least some
% well-understood assumption horizontally?
% 2) why is the cone density different in the left and right eyes in
% ISETBio? 
% 3) and a note: title of plot is no longer accurate, has not been updated
% since both species were plotted


%choose size of screen
screen = -10:.1:10;

D = 10; %cm distance from eyes to screen

distanceFalloff = 1; %how much does information content fall off with distance to object

%TS
phiTS = 20; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
E_dTS = 0; % cm, distance between eyes

%H
phiH = 0;
E_dH = 0;

%optics parameters:
tsSigma = 3; % max cone density is 36k at center, but how long before it drops to 12k?

% initate result vectors
leftTSEcc = [];
leftTSConeDensity = [];
leftTSDistance = [];
leftTSI = [];
rightTSEcc = [];
rightTSConeDensity = [];
rightTSDistance = [];
rightTSI = [];
leftHumanEcc = [];
leftHumanConeDensity = [];
leftHumanDistance = [];
leftHumanI = [];
rightHumanEcc = [];
rightHumanConeDensity = [];
rightHumanDistance = [];
rightHumanI = [];
sumHumanI = [];
sumTSI = [];


for i = 1:length(screen)

%TS    
leftTSEcc(i) = atand((screen(i)-(E_dTS/2))/D)-phiTS;
leftTSConeDensity(i) = normpdf(leftTSEcc(i),0,tsSigma)*(24000/.4)+12000;
leftTSDistance(i) = D/cosd(phiTS+leftTSEcc(i));
leftTSI(i) = leftTSConeDensity(i)/(leftTSDistance(i)^1);

rightTSEcc(i) = atand((-screen(i)-(E_dTS/2))/D)-phiTS; 
rightTSConeDensity(i) = normpdf(rightTSEcc(i),0,tsSigma)*(24000/.4)+12000;
rightTSDistance(i) = D/cosd(phiTS+rightTSEcc(i));
rightTSI(i) = rightTSConeDensity(i)/(rightTSDistance(i)^1);

%human
leftHumanEcc(i) = atand((screen(i)-(E_dH/2))/D)-phiH;
leftHumanConeDensity(i) = coneDensityReadData('eccentricity', abs((leftHumanEcc(i))*(3/1000)), 'whichEye', 'left');
if(isnan(leftHumanConeDensity(i)))
    leftHumanConeDensity(i) = 4.8145e+03;
end
leftHumanDistance(i) = D/cosd(phiH+leftHumanEcc(i));
leftHumanI(i) = leftHumanConeDensity(i)/(leftHumanDistance(i)^1);

rightHumanEcc(i) = atand((screen(i)-(E_dH/2))/D)-phiH;
rightHumanConeDensity(i) = coneDensityReadData('eccentricity', abs((rightHumanEcc(i))*(3/1000)), 'whichEye', 'right');
if(isnan(rightHumanConeDensity(i)))
    rightHumanConeDensity(i) = 3.2540e+03;
end
rightHumanDistance(i) = D/cosd(phiH+rightHumanEcc(i));
rightHumanI(i) = rightHumanConeDensity(i)/(rightHumanDistance(i)^1);

%%
sumHumanI(i) = leftHumanI(i) + rightHumanI(i);
sumTSI(i) = leftTSI(i) + rightTSI(i);

end

%hold on
hold on
plot(screen,sumTSI)
plot(screen,sumHumanI)
xlabel('Horizontal Meridian (cm)')
ylabel('Cone Density Combined Across Eyes')
title(sprintf('Information Effiency Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
hold off
%set(gca, 'XLabel', 'Horizontal Meridian (cm)','YLabel','Cone Density Combined Across Eyes')
