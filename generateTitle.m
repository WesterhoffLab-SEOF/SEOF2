
function title = generateTitle(SystemParam)
    title = sprintf('Irradiance along the Fiber Length\n %d Fiber(s), %d cm long, %d [μm] OD , %d mm coupling distance', ...
        SystemParam.numFibers, SystemParam.xLen/1e4, SystemParam.fiberRadius*2, SystemParam.ledDistance);
end