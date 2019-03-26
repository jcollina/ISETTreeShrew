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

screenSize = 15;

D = 1; %cm distance from eyes to screen

distanceFalloff = 1; %how much does information content fall off with distance to object

human = 0;
treeShrew = 1;
mouse = 0;

if human
    
    phiH = 0;
    E_dH = 0;
    
    leftHumanEcc = zeros(size(screen));
    leftHumanConeDensity = zeros(size(screen));
    leftHumanDistance = zeros(size(screen));
    leftHumanI = zeros(size(screen));
    rightHumanEcc = zeros(size(screen));
    rightHumanConeDensity = zeros(size(screen));
    rightHumanDistance = zeros(size(screen));
    rightHumanI = zeros(size(screen));
    
    sumConeDensityHuman = zeros(size(screen));
    sumHumanI = zeros(size(screen));
end

if treeShrew
    phiTS = 0; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
    E_dTS = 0; % cm, distance between eyes
    
    leftTSEcc = zeros(size(screen));
    leftTSConeDensity = zeros(size(screen));
    leftTSDistance = zeros(size(screen));
    leftTSI = zeros(size(screen));
    rightTSEcc = zeros(size(screen));
    rightTSConeDensity = zeros(size(screen));
    rightTSDistance = zeros(size(screen));
    rightTSI = zeros(size(screen));
    
    sumConeDensityTS = zeros(size(screen));
    sumTSI = zeros(size(screen));
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

tsSigma = 3; %degrees... max cone density is 36k at center, but how long before it drops to 12k?

axialLengthHuman = 23; %mm
axialLengthTS = 8; %mm
axialLengthMouse = 3.62; %mm


totalHumanFOV = 120;
totalTSFOV = 120;
totalMouseFOV = 120;


maxEccTS = totalTSFOV/2;
maxEccHuman = totalHumanFOV/2;
maxEccMouse = totalMouseFOV/2;
%%
screen = cell(10);
for ii = 1:size(screen,1)
    for jj = 1:size(screen,2)
        c = [(size(screen,1)-ii)-(size(screen,1)/2-0.5),(jj-1)-(size(screen,2)/2-0.5)];
        screen{ii,jj}=c;
    end
end
%%
for x = 1:size(screen,1) % x
   for y = 1:size(screen,1) % y 


    coords = screen{x,y};
    screenX = coords(1)
    screenY = coords(2)
    
    %if treeShrew
    if screenY >= 0
        angle = atand(screenY/screenX);
    else
        angle = 180+atand(screenY/screenX);
    end
        actualDistance = sqrt(sqrt(D^2+screenX^2)+screenY^2);
        observedDistance = sqrt(sqrt(D^2+(screenX+D*cos(phiTS))^2)+screenY^2);
        TSEcc = atand(screenY/observedDistance);

        if abs(TSEcc(x,y)) > maxEccTS
            leftTSConeDensity(x,y) = 0;
        else
            leftTSConeDensity(x,y) = normpdf(leftTSEcc(x,y),0,tsSigma)*(24000/.4)+12000;
        end
        
        %leftTSDistance(x,y) = D/cosd(phiTS+leftTSEcc(x,y));
        
        
        leftTSI(x,y) = leftTSConeDensity(x,y)/(leftTSDistance(x,y)^1);
        
        
        rightTSEcc(x,y) = atand((-screenX-(E_dTS/2))/D)+phiTS;
        if abs(rightTSEcc(x,y)) > maxEccTS
            rightTSConeDensity(x,y) = 0;
        else
            rightTSConeDensity(x,y) = normpdf(rightTSEcc(x,y),0,tsSigma)*(24000/.4)+12000;
        end
        
        rightTSDistance(x,y) = D/cosd(phiTS+rightTSEcc(x,y));
        
        rightTSI(x,y) = rightTSConeDensity(x,y)/(rightTSDistance(x,y)^1);
        
        sumConeDensityTS(x,y) = leftTSConeDensity(x,y) + rightTSConeDensity(x,y);
        
        sumTSI(x,y) = (leftTSI(x,y) + rightTSI(x,y))*(axialLengthTS/2)^2;
    end
    
   end
%end
%{    
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
%}
%hold on
figure(1)
hold on
if treeShrew
    plot(screen,sumConeDensityTS)
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
