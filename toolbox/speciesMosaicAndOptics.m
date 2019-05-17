%%
function theOI = speciesOI(species,varargin)
disp(varargin)
    switch species
        case "human"
            theOI = oiCreate(varargin{:});
        case "treeshrew"
            theOI = oiTreeShrewCreate(varargin{:});
        otherwise
            error("Species should be treeshrew or human... for now")
    end   
end
%% ewe
function theMosaic = speciesMosaic(species,varargin)
    switch species
        case "human"
            theMosaic = coneMosaicHex(7,varargin{:});
        case "treeshrew"
            theMosaic = coneMosaicTreeShrewCreate(75,varargin{:});
        otherwise
            error("Species should be treeshrew or human... for now")
    end
end
%{


            theMosaic = coneMosaicTreeShrewCreate(75, ...
                'fovDegs', sizeDegs, ...
                'customLambda', cone_spacing);
                            
            theMosaic = coneMosaicHex(7, ...
                'eccBasedConeDensity', true, ...
                'maxGridAdjustmentIterations',200, ...
                'fovDegs', sizeDegs);
%}