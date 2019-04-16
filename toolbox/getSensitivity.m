function [isetSensitivity,geomSensitivity] = getSensitivity(pupilDiameterMM,focalLengthMM,innerSegmentDiameterUM,testScene,spatialLMSdensities,fovDegs,whichConeType);        

% Create optical image object for current parameters.
    tOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM, 'focalLengthMM', ...
        focalLengthMM);
    
    % Create the mosaic
    tMosaic = coneMosaicTreeShrewCreate(tOI.optics.micronsPerDegree, ...
        'spatialDensity', spatialLMSdensities, ...
        'customInnerSegmentDiameter', innerSegmentDiameterUM, ...
        'integrationTime', 5/1000, ...
        'fovDegs', fovDegs);
    
    % Compute the retinal image
    tOI = oiCompute(tOI, testScene);
    
    % Get and store the retinal irradiance
    meanRetinalIlluminance = oiGet(tOI, 'mean illuminance');
    
    % Compute the mosaic responses (for more replicable responses, use more trials)
    nTrialsNum = 1;
    emPath = zeros(nTrialsNum, 1, 2);
    
    % Compute *treeshrew* mosaic excitation responses to treeshrew optical image
    tMosaicExcitation = tMosaic.compute(tOI, 'emPath', emPath);
    
    % Find mean excitations
    isetSensitivity = ...
        meanResponseToOpticalImage(tMosaic, tMosaicExcitation, whichConeType);
    
    % Calculate estimated sensitivity according to Animal Eyes
    geomSensitivity = 0.62 * (pupilDiameterMM^2 * innerSegmentDiameterUM^2)/ ...
        (focalLengthMM^2);
end

function meanResponse = meanResponseToOpticalImage(coneMosaic, coneMosaicResponse, ...
    targetConeType)
nTrialsNum = size(coneMosaicResponse,1);
coneMosaicResponse  = reshape(coneMosaicResponse, [nTrialsNum numel(coneMosaic.pattern)]);
idx = find(coneMosaic.pattern == targetConeType);
meanResponse = mean(mean(coneMosaicResponse(:,idx)));
end