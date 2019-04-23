b = 11:2:15;
timeS=[];
mosTs=[];
for s = 1:length(b)
    
    sizeDegs = b(s);
    
    %
    nTrialsNum = 100;
    tm = tic;
    theMosaic = coneMosaicTreeShrewCreate(75, ...%theOI.optics.micronsPerDegree, ...
        'fovDegs', sizeDegs, ...        % match mosaic width to stimulus size
        'integrationTimeSeconds', 10/1000);
    mosTs(s) = toc(tm);
    %
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);
    
    % parameter struct for a low spatial frequency Gabor stimulus
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', 1, ...  % changing cycles/deg
        'orientationDegs', 0, ...               % 0 degrees
        'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
        'sizeDegs', sizeDegs, ...               % 14 x 14 size
        'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
        'contrast', 0.1,...                       % 0 Michelson contrast
        'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
        'pixelsAlongWidthDim', [], ...          % pixels- width dimension
        'pixelsAlongHeightDim', [] ...          % pixel- height dimension
        );
    
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);
    
    stimParams.contrast = 0.0;
    
    nullScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);
    [~,time] = getSVMAcc(theMosaic, testScene, nullScene, nTrialsNum);
    timeS(s) = time;
    
end
mosTs
timeS