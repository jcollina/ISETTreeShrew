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

compute = 1;
plotting = 1;

human = 1;
treeShrew = 1;
mouse = 0;

%choose size of screen

screenSize = 100/1000; %meters

screenDistance = 100/1000; %meters distance from eyes to screen

distanceFalloff = .25; %how fast does information content fall off with distance to object?

%optics parameters:

tsSigma = 1; %degrees... max cone density is 36k at center, but how long before it drops to 12k?

axialLengthH = 23/1000; %m
axialLengthTS = 8/1000; %m
axialLengthM = 3.62/1000; %m


totalHumanFOV = 120; %deg
totalTSFOV = 120;
totalMouseFOV = 120;

if compute
    
    dimScreen = screenSize * 1000 + 1;
    
    D = screenDistance * 1000;
    
    maxEccTS = totalTSFOV/2;
    maxEccH = totalHumanFOV/2;
    maxEccM = totalMouseFOV/2;
    
    if human
        phiH = 0; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
        E_d_H = 0/1000; % m, distance between eyes
        
        L_CD_H = zeros(dimScreen);
        R_CD_H = zeros(dimScreen);
        
        L_I_H = zeros(dimScreen);
        R_I_H = zeros(dimScreen);
        
        sum_CD_H = zeros(dimScreen);
        sum_I_H = zeros(dimScreen);
    end
    
    if treeShrew
        phiTS = 25; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
        E_d_TS = 0/1000; % m, distance between eyes
        
        L_CD_TS = zeros(dimScreen);
        R_CD_TS = zeros(dimScreen);
        
        L_I_TS = zeros(dimScreen);
        R_I_TS = zeros(dimScreen);
        
        sum_CD_TS = zeros(dimScreen);
        sum_I_TS = zeros(dimScreen);
    end
    
    if mouse
        phiM = 0; %degrees eye is rotated (0 degrees: fovea pointed at screen, 90 degrees: fovea pointed adjacent)
        E_d_M = 20/1000; % m, distance between eyes
        
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
            screenTemp{ii,jj}=c/1000;
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
                    L_CD_TS = 0;
                else
                    L_CD_TS = normpdf(ecc(1),0,tsSigma)*(24000/.4)+12000;
                end
                
                L_I_TS(ii,jj) = L_CD_TS/(dist(1)^1);
                
                if abs(ecc(2)) > maxEccTS
                    R_CD_TS = 0;
                else
                    R_CD_TS = normpdf(ecc(2),0,tsSigma)*(24000/.4)+12000;
                end
                
                R_I_TS(ii,jj) = R_CD_TS/(dist(2)^1);
                
                sum_CD_TS(ii,jj) = L_CD_TS + R_CD_TS;
                
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
                L_I_H(ii,jj) = L_CD_H/(dist(1)^1);
                
                if abs(ecc(2)) > maxEccH
                    R_CD_H = 0;
                else
                    R_CD_H = coneDensityReadData('eccentricity', abs(ecc(2))*(3/1000),'angle',angle(2), 'whichEye', 'right');
                    if(isnan(R_CD_H))
                        R_CD_H = 3.2540e+03;
                    end
                end
                
                R_I_H(ii,jj) = R_CD_H/(dist(2)^1);
                
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
                
                L_I_M(ii,jj) = L_CD_M/(dist(1)^1);
                
                if abs(ecc(2)) > maxEccM
                    R_CD_M = 0;
                else
                    R_CD_M = 12000;%normpdf(ecc(2),0,tsSigma)*(24000/.4)+12000;
                end
                
                R_I_M(ii,jj) = R_CD_M/(dist(2)^1);
                
                sum_CD_M(ii,jj) = L_CD_M + R_CD_M;
                
                sum_I_M(ii,jj) = (L_I_M(ii,jj) + R_I_M(ii,jj))*(axialLengthM/2)^2;
            end
        end
    end
end

if(plotting)
    
    %hold on
    s = screen; %to get in mm for plotting
    fun = @(x) 1000 * x(1);
    X = cellfun(fun,s);
    fun = @(x) 1000 * x(2);
    Y = cellfun(fun,s);
    
    figure(1)
    
    if human
        c = getColMat(dimScreen,[0,1,0]);
        surf(X,Y,log10(sum_CD_H),'FaceColor','g', 'FaceAlpha',0.3, 'EdgeColor','none')
        hold on
    end
    
    if treeShrew
        c = getColMat(dimScreen,[1,0,0]);
        surf(X,Y,log10(sum_CD_TS),'FaceColor','r', 'FaceAlpha',0.3)%, 'EdgeColor','none')
        hold on
    end
    
    if mouse
        c = getColMat(dimScreen,[0,0,1]);
        surf(X,Y,sum_I_M,'FaceColor','b', 'FaceAlpha',0.3)%, 'EdgeColor','none')
        hold on
    end
    xlabel('Displacement from Horizontal Meridian (mm)')
    ylabel('Displacement from Vertical Meridian (mm)')
    zlabel('log_{10} Cone Density (cones/mm^2)')
    %title(sprintf('Cone Density Due to Eye Angle \n D=%.0f cm, \\phi=%.0f degrees, E_d = %.0f cm', D,phi,E_dTS))
    hold off
    grid on
    
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