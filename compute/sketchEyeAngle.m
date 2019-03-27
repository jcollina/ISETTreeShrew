% questions to ask:
% 1) what's known about the specific drop-off in cone density? It
% definitely doesn't seem spherically symmetric, but is there at least some
% well-understood assumption horizontally?
% 2) why is the cone density different in the left and right eyes in
% ISETBio?
% 3) and a note: title of plot is no longer accurate, has not been updated
% since both species were plotted

%note: used cone density from humans from ISETBio, made guesses about TS
%using knowledge that TS density ranges from 12,000 to 36,000. Also made
%guesses about mouse- no fovea, and average of 12,000 cones/mm^2, so just
%had that be the cone density everywhere that was in the field of vision

%currently have max FOV for all eyes being 120 total degrees (so 60 each
%side)

%NEXT STEP: bring in size of eyes: break up "visual field" into degrees-
%it's not just about the density, but about the number of cones that the
%eye has on that area.

%surface area of sphere varies with radius as r^2, so for approx, just
%multiply cone density by ((axial length)/2)

%if all we care about is cone density, then we can see clear areas in the visual field where
%the tree shrew AND mouse eyes have more cones/mm^2 than humans. But humans
%have much bigger eyes, so even in those areas, there are more total cones
%in the human eye for each segment of visual area.


%choose size of screen
screen = -100:.1:100;

D = 10; %cm distance from eyes to screen

distanceFalloff = 1; %how much does information content fall off with distance to object

human = 0;
treeShrew = 1;
mouse = 0;

if human
    phiH = 0;
    E_dH = 0;
    
    leftHumanEcc = zeros(1,length(screen));
    leftHumanConeDensity = zeros(1,length(screen));
    leftHumanDistance = zeros(1,length(screen));
    leftHumanI = zeros(1,length(screen));
    rightHumanEcc = zeros(1,length(screen));
    rightHumanConeDensity = zeros(1,length(screen));
    rightHumanDistance = zeros(1,length(screen));
    rightHumanI = zeros(1,length(screen));
    
    sumConeDensityHuman = zeros(1,length(screen));
    sumHumanI = zeros(1,length(screen));
end

if treeShrew
    phiTS = 20; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
    E_dTS = 0; % cm, distance between eyes
    
    leftTSEcc = zeros(1,length(screen));
    leftTSConeDensity = zeros(1,length(screen));
    leftTSDistance = zeros(1,length(screen));
    leftTSI = zeros(1,length(screen));
    rightTSEcc = zeros(1,length(screen));
    rightTSConeDensity = zeros(1,length(screen));
    rightTSDistance = zeros(1,length(screen));
    rightTSI = zeros(1,length(screen));
    
    sumConeDensityTS = zeros(1,length(screen));
    sumTSI = zeros(1,length(screen));
end

if mouse
    phiM = 20; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
    E_dM = 0; % cm, distance between eyes
    
    leftMouseEcc = zeros(1,length(screen));
    leftMouseConeDensity = zeros(1,length(screen));
    leftMouseDistance = zeros(1,length(screen));
    leftMouseI = zeros(1,length(screen));
    rightMouseEcc = zeros(1,length(screen));
    rightMouseConeDensity = zeros(1,length(screen));
    rightMouseDistance = zeros(1,length(screen));
    rightMouseI = zeros(1,length(screen));
    
    sumConeDensityMouse = zeros(1,length(screen));
    sumMouseI = zeros(1,length(screen));
end

%optics parameters:

tsSigma = 10; %degrees... max cone density is 36k at center, but how long before it drops to 12k?

axialLengthHuman = 23; %mm
axialLengthTS = 8; %mm
axialLengthMouse = 3.62; %mm


totalHumanFOV = 120;
totalTSFOV = 1120;
totalMouseFOV = 120;


maxEccTS = totalTSFOV/2;
maxEccHuman = totalHumanFOV/2;
maxEccMouse = totalMouseFOV/2;

for i = 1:length(screen)
    if treeShrew
        %TS
        leftTSEcc(i) = atand((screen(i)-(E_dTS/2))/D)+phiTS;
        
        if abs(leftTSEcc(i)) > maxEccTS
            leftTSConeDensity(i) = 0;
        else
            leftTSConeDensity(i) = normpdf(leftTSEcc(i),0,tsSigma)*(24000/.4)+12000;
        end
        
        leftTSDistance(i) = D/cosd(phiTS+leftTSEcc(i));
        
        leftTSI(i) = leftTSConeDensity(i)/(leftTSDistance(i)^1);
        
        rightTSEcc(i) = atand((-screen(i)-(E_dTS/2))/D)+phiTS;
        if abs(rightTSEcc(i)) > maxEccTS
            rightTSConeDensity(i) = 0;
        else
            rightTSConeDensity(i) = normpdf(rightTSEcc(i),0,tsSigma)*(24000/.4)+12000;
        end
        
        rightTSDistance(i) = D/cosd(phiTS+rightTSEcc(i));
        
        rightTSI(i) = rightTSConeDensity(i)/(rightTSDistance(i)^1);
        
        sumConeDensityTS(i) = leftTSConeDensity(i) + rightTSConeDensity(i);
        
        sumTSI(i) = (leftTSI(i) + rightTSI(i))*(axialLengthTS/2)^2;
    end
    
    if human
        %human
        leftHumanEcc(i) = atand((screen(i)-(E_dH/2))/D)+phiH;
        
        if abs(leftHumanEcc(i)) > maxEccHuman
            leftHumanConeDensity(i) = 0;
        else
            leftHumanConeDensity(i) = coneDensityReadData('eccentricity', abs((leftHumanEcc(i))*(3/1000)), 'whichEye', 'left');
        end
        
        if(isnan(leftHumanConeDensity(i)))
            leftHumanConeDensity(i) = 4.8145e+03; %would love a better way of doing this- maybe read in all data first, then get lowest, and then us that?
        end
        
        leftHumanDistance(i) = D/cosd(phiH+leftHumanEcc(i));
        
        leftHumanI(i) = leftHumanConeDensity(i)/(leftHumanDistance(i)^1);
        
        
        rightHumanEcc(i) = atand((screen(i)-(E_dH/2))/D)+phiH;
        if abs(rightHumanEcc(i)) > maxEccHuman
            rightHumanConeDensity(i) = 0;
        else
            rightHumanConeDensity(i) = coneDensityReadData('eccentricity', abs((rightHumanEcc(i))*(3/1000)), 'whichEye', 'right');
        end
        if(isnan(rightHumanConeDensity(i)))
            rightHumanConeDensity(i) = 3.2540e+03;
        end
        
        rightHumanDistance(i) = D/cosd(phiH+rightHumanEcc(i));
        
        
        rightHumanI(i) = rightHumanConeDensity(i)/(rightHumanDistance(i)^1);
        
        sumConeDensityHuman(i) = leftHumanConeDensity(i) + rightHumanConeDensity(i);
        sumHumanI(i) = (leftHumanI(i) + rightHumanI(i))*(axialLengthHuman/2)^2;
    end
    
    
    if mouse
        
        leftMouseEcc(i) = atand((screen(i)-(E_dM/2))/D)+phiM;
        
        if abs(leftMouseEcc(i)) > maxEccMouse
            leftMouseConeDensity(i) = 0;
        else
            leftMouseConeDensity(i) = 12000;
        end
        
        leftMouseDistance(i) = D/cosd(phiM+leftMouseEcc(i));
        
        leftMouseI(i) = leftMouseConeDensity(i)/(leftMouseDistance(i)^1);
        
        rightMouseEcc(i) = atand((-screen(i)-(E_dM/2))/D)+phiM;
        
        rightMouseDistance(i) = D/cosd(phiM+rightMouseEcc(i));
        
        if abs(rightMouseEcc(i)) > maxEccMouse
            rightMouseConeDensity(i) = 0;
        else
            rightMouseConeDensity(i) = 12000;
        end
        
        rightMouseI(i) = rightMouseConeDensity(i)/(rightMouseDistance(i)^1);
        
        sumConeDensityMouse(i) = leftMouseConeDensity(i) + rightMouseConeDensity(i);
        
        sumMouseI(i) = (leftMouseI(i) + rightMouseI(i))*(axialLengthMouse/2)^2;
    end

end

%hold on
figure(1)
hold on
if treeShrew
    plot(screen,leftTSConeDensity)
end
if human
    plot(screen,sumConeDensityHuman)
end
if mouse
    plot(screen,sumConeDensityMouse)
end
xlabel('Horizontal Meridian (cm)')
ylabel('Cone Density Combined Across Eyes')
%title(sprintf('Information Effiency Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
hold off
%set(gca, 'XLabel', 'Horizontal Meridian (cm)','YLabel','Cone Density Combined Across Eyes')

figure(2)
hold on
if treeShrew
    plot(screen,sumTSI)
end
if human
    plot(screen,sumHumanI)
end
if mouse
    plot(screen,sumMouseI)
end
xlabel('Horizontal Meridian (cm)')
ylabel('Information Combined Across Eyes')
%title(sprintf('Information Effiency Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
hold off
