function theOI = speciesOI(species,varargin)
class(varargin)
    switch species
        case "human"
            theOI = oiCreate(varargin{:});
        case "treeshrew"
            theOI = oiTreeShrewCreate(varargin{:});
        otherwise
            error("Species should be treeshrew or human... for now")
    end   
end