%% t_eyeGeometry
% Illustrate effect of various geometric properties of the eye on
% information obtained from an image
%

%
% See also:
%

% History:
%   03/26/19 jsc  Wrote initial version.

human = 1;
treeShrew = 1;
mouse = 1;

compute = 1;

% do I want:
% just the plot of cone density for each spot on a screen
coneDensityPlot = 1;

%here, I take into account information falloff due to distance, and information
%increase due to eye size
eyeInformationPlot = 1;

%choose size and distance of screen

screenSize = 1000/1000; %meters

screenDistance = 500/1000; %meters distance from eyes to screen

%optics parameters:

distanceFalloff = .25; %how fast does information content fall off with distance to object?

axialLengthH = 23/1000; %meters
phiH = 0; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
E_d_H = 0/1000; % meters, distance between eyes
totalHumanFOV = 120; %deg

axialLengthTS = 8/1000;
phiTS = 30;
E_d_TS = 0/1000;
totalTSFOV = 120;
tsSigma = 10; %degrees... max cone density is 36k at center, but how long before it drops to 12k?

axialLengthM = 3.62/1000;
phiM = 30;
E_d_M = 20/1000;
totalMouseFOV = 120;


if compute
    
    dimScreen = 101;%screenSize * 1000 + 1;
    
    D = screenDistance * 1000;
    
    maxEccTS = totalTSFOV/2;
    maxEccH = totalHumanFOV/2;
    maxEccM = totalMouseFOV/2;
    
    maxP = normpdf(0,0,tsSigma);
    
    if human
        
        L_CD_H = zeros(dimScreen);
        R_CD_H = zeros(dimScreen);
        
        L_I_H = zeros(dimScreen);
        R_I_H = zeros(dimScreen);
        
        sum_CD_H = zeros(dimScreen);
        sum_I_H = zeros(dimScreen);
    end
    
    if treeShrew
        
        L_CD_TS = zeros(dimScreen);
        R_CD_TS = zeros(dimScreen);
        
        L_I_TS = zeros(dimScreen);
        R_I_TS = zeros(dimScreen);
        
        sum_CD_TS = zeros(dimScreen);
        sum_I_TS = zeros(dimScreen);
    end
    
    if mouse
        
        
        L_CD_M = zeros(dimScreen);
        R_CD_M = zeros(dimScreen);
        
        L_I_M = zeros(dimScreen);
        R_I_M = zeros(dimScreen);
        
        sum_CD_M = zeros(dimScreen);
        sum_I_M = zeros(dimScreen);
    end
    
    
    %
    
    screenTemp = cell(dimScreen);
    
    for ii = 1:size(screenTemp,1)
        for jj = 1:size(screenTemp,2)
            c = [(jj-1)-(size(screenTemp,2)/2-0.5),(size(screenTemp,1)-ii)-(size(screenTemp,1)/2-0.5)];
            screenTemp{ii,jj} = (c/50) * screenSize/2;
        end
    end
    
    screen = screenTemp;
    
    %
    for ii = 1:size(screen,1) % x
        for jj = 1:size(screen,1) % z
            
            %%
            coords = screen{ii,jj};
            
            if treeShrew

                [dist,obsDist,ecc,angle] = getOpticsGeometry(coords,screenDistance,phiTS,E_d_TS);
                angles(ii,jj) = angle(1);
                eccs(ii,jj) = ecc(1);
                
                if abs(ecc(1)) > maxEccTS
                    L_CD_TS(ii,jj) = 0;
                else
                    L_CD_TS(ii,jj) = normpdf(ecc(1),0,tsSigma)*(24000/maxP)+12000;
                end
                
                L_I_TS(ii,jj) = L_CD_TS(ii,jj)/(dist(1)^distanceFalloff);
                
                if abs(ecc(2)) > maxEccTS
                    R_CD_TS = 0;
                else
                    R_CD_TS = normpdf(ecc(2),0,tsSigma)*(24000/maxP)+12000;
                end
                
                R_I_TS(ii,jj) = R_CD_TS/(dist(2)^distanceFalloff);
                
                sum_CD_TS(ii,jj) = L_CD_TS(ii,jj) + R_CD_TS;
                
                sum_I_TS(ii,jj) = (L_I_TS(ii,jj) + R_I_TS(ii,jj))*(axialLengthTS/2)^2;
            end
            
            
            
            %end
            
            if human
                %human
                
                [dist,obsDist,ecc,angle] = getOpticsGeometry(coords,screenDistance,phiH,E_d_H);
                angles(ii,jj) = angle(1);
                eccs(ii,jj) = ecc(1);
                if abs(ecc(1)) > maxEccH
                    L_CD_H = 0;
                else
                    L_CD_H = coneDensityReadData('eccentricity', abs(ecc(1))*(3/1000),'angle',angle(1), 'whichEye', 'left');
                    if(isnan(L_CD_H))
                        L_CD_H = 4.8145e+03; %would love a better way of doing this- maybe read in all data first, then get lowest, and then us that?
                    end
                end
                
                ls(ii,jj)=L_CD_H;
                dists(ii,jj)=dist(1);
                L_I_H(ii,jj) = L_CD_H/(dist(1)^distanceFalloff);
                
                if abs(ecc(2)) > maxEccH
                    R_CD_H = 0;
                else
                    R_CD_H = coneDensityReadData('eccentricity', abs(ecc(2))*(3/1000),'angle',angle(2), 'whichEye', 'right');
                    if(isnan(R_CD_H))
                        R_CD_H = 3.2540e+03;
                    end
                end
                
                R_I_H(ii,jj) = R_CD_H/(dist(2)^distanceFalloff);
                
                sum_CD_H(ii,jj) = L_CD_H + R_CD_H;
                
                sum_I_H(ii,jj) = (L_I_H(ii,jj) + R_I_H(ii,jj))*(axialLengthH/2)^2;
                
            end
            
            
            if mouse
                
                [dist,obsDist,ecc,angle] = getOpticsGeometry(coords,screenDistance,phiM,E_d_M);
                
                if abs(ecc(1)) > maxEccM
                    L_CD_M = 0;
                else
                    L_CD_M = 12000;%normpdf(ecc(1),0,tsSigma)*(24000/.4)+12000;
                end
                
                L_I_M(ii,jj) = L_CD_M/(dist(1)^distanceFalloff);
                
                if abs(ecc(2)) > maxEccM
                    R_CD_M = 0;
                else
                    R_CD_M = 12000;%normpdf(ecc(2),0,tsSigma)*(24000/.4)+12000;
                end
                
                R_I_M(ii,jj) = R_CD_M/(dist(2)^distanceFalloff);
                
                sum_CD_M(ii,jj) = L_CD_M + R_CD_M;
                
                sum_I_M(ii,jj) = (L_I_M(ii,jj) + R_I_M(ii,jj))*(axialLengthM/2)^2;
            end
        end
    end
end

if(coneDensityPlot)
    
    %hold on
    s = screen; %to get in mm for plotting
    fun = @(x) 1000 * x(1);
    X = cellfun(fun,s);
    fun = @(x) 1000 * x(2);
    Y = cellfun(fun,s);
    
    dimMax = max(X(:));
    
    figure(1)
    
    if human
        Z = log10(sum_CD_H);
        surf(X,Y,Z,'FaceColor','r', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'r-');
        plot3(X, 0*Y + dimMax, Z,'r.-');
        plot3(0*X + dimMax, Y, Z,'r.-'); 
    end
    if treeShrew
        Z = log10(sum_CD_TS);
        surf(X,Y,Z,'FaceColor','g', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'g-');
        plot3(X, 0*Y + dimMax, Z,'g.-');
        plot3(0*X + dimMax, Y, Z,'g.-'); 
    end
    if mouse
        Z = log10(sum_CD_M);
        surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'b-');
        plot3(X, 0*Y + dimMax, Z,'b.-');
        plot3(0*X + dimMax, Y, Z,'b.-'); 
    end
    xlabel('Displacement from Horizontal Meridian (mm)')
    ylabel('Displacement from Vertical Meridian (mm)')
    zlabel('log_{10} Cone Density (cones/mm^2)')
    %title(sprintf('Cone Density Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
    hold off
    grid on
   
end

if eyeInformationPlot
    
    s = screen; %to get in mm for plotting
    fun = @(x) 1000 * x(1);
    X = cellfun(fun,s);
    fun = @(x) 1000 * x(2);
    Y = cellfun(fun,s);
    
    figure(2)
    
    if human
        Z = log10(sum_I_H);
        surf(X,Y,Z,'FaceColor','r', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'r-');
        plot3(X, 0*Y + dimMax, Z,'r.-');
        plot3(0*X + dimMax, Y, Z,'r.-'); 
    end
    if treeShrew
        Z = log10(sum_I_TS);
        surf(X,Y,Z,'FaceColor','g', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'g-');
        plot3(X, 0*Y + dimMax, Z,'g.-');
        plot3(0*X + dimMax, Y, Z,'g.-'); 
    end
    if mouse
        Z = log10(sum_I_M);
        surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
        %plot3(X, Y, 0*Z + dimMax,'b-');
        plot3(X, 0*Y + dimMax, Z,'b.-');
        plot3(0*X + dimMax, Y, Z,'b.-'); 
    end
    xlabel('Displacement from Horizontal Meridian (mm)')
    ylabel('Displacement from Vertical Meridian (mm)')
    zlabel('log_{10} Information (~Cones)')
    %title(sprintf('Cone Density Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
    hold off
    grid on
end

%% Functions

function [distance,observedDistance,eccentricity,angle] = getOpticsGeometry(coords,D,phi,E_d)
ScreenZ = coords(2);
screenPosVector = [coords(1) 0 coords(2)];
b = [-1 1]; % left, right
for i = 1:length(b)
    side = b(i);
    eyePosVector = [side * E_d/2 D 0];
    screenX = screenPosVector(1)-eyePosVector(1);
    observedScreenX = screenX - side * D * tand(phi);
    observedScreenPosVector = [observedScreenX 0 coords(2)];
    if observedScreenX >= 0
        angle(i) = atand(ScreenZ/(observedScreenX));
    else
        angle(i) = 180+atand(ScreenZ/observedScreenX);
    end
    if isnan(angle(i))
        angle(i) = 1e-4;
    end
    distance(i) = sqrt(sum((eyePosVector - screenPosVector).^2));
    observedDistance(i) = sqrt(sum((eyePosVector - observedScreenPosVector).^2));
    eccentricity(i) = atand((sqrt(observedScreenX^2+ScreenZ^2))/distance(i));
end
end

function c = getColMat(screenSize,rgb)
if(prod(rgb)>1)
    rgb = rgb./255;
end
c = zeros(screenSize,screenSize,3);
c(:,:,1) = rgb(1);
c(:,:,2) = rgb(2);
c(:,:,3) = rgb(3);
end