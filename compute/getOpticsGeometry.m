function [distance,observedDistance,eccentricity,angle] = getOpticsGeometry(coords,D,phi,E_d)

ScreenZ = coords(2);

screenPosVector = [coords(1) 0 coords(2)];

%left
b = [-1 1];
for i = 1:length(b)

    
    side = b(i); %left/right
    eyePosVector = [side * E_d/2 D 0];
    
    screenX = screenPosVector(1)-eyePosVector(1);
    
    observedScreenPosVector = [coords(1) - side * D * tand(phi) 0 coords(2)];
    
    if ScreenZ >= 0
        angle(i) = atand(ScreenZ/(screenX));
    else
        angle(i) = 180+atand(ScreenZ/screenX);
    end
    
    distance(i) = sqrt(sum((eyePosVector - screenPosVector).^2));
    observedDistance(i) = sqrt(sum((eyePosVector - observedScreenPosVector).^2));
    eccentricity(i) = atand(ScreenZ/observedDistance(i));
end

end