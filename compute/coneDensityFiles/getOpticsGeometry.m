function [distance,observedDistance,eccentricity,angle] = getOpticsGeometry(coords,D,phi,E_d)

ScreenZ = coords(2);

screenPosVector = [coords(1) 0 coords(2)];

%left
b = [-1 1];
for i = 1:length(b)

    
    side = b(i); %left/right
    eyePosVector = [side * E_d/2 D 0];
    
    screenX = screenPosVector(1)-eyePosVector(1);
    observedScreenX(i) = screenX - side * D * tand(phi);
    observedScreenPosVector = [observedScreenX(i) 0 coords(2)];
    
    if observedScreenX >= 0
        angle(i) = atand(ScreenZ/(observedScreenX(i)));
    else
        angle(i) = 180+atand(ScreenZ/observedScreenX(i));
    end
    
    if isnan(angle(i))
        angle(i) = 1e-4;
    end
    
    distance(i) = sqrt(sum((eyePosVector - screenPosVector).^2));
    observedDistance(i) = sqrt(sum((eyePosVector - observedScreenPosVector).^2));
    
    eccentricity(i) = atand((sqrt(observedScreenX(i)^2+ScreenZ^2))/distance(i));
end
%%
end