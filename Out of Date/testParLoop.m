presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 1, ...  % changing cycles/deg
    'contrast',.5, ...
    'orientationDegs', 0, ...               % 0 degrees
    'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
    'sizeDegs', 1, ...               % size
    'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
    'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
    'pixelsAlongWidthDim', [], ...          % pixels- width dimension
    'pixelsAlongHeightDim', [] ...          % pixel- height dimension
);

parfor i = 1:4
    %a = load('luminosity.mat');
    %disp(a)
    isetbioRootPath()
    for a = 1:100
        b = a * 10
    end
    nullScene = generateGaborScene(...
        "stimParams", stimParams,...
        "presentationDisplay", presentationDisplay);
    
end