classdef travel_storage
    properties
        absorbi;%sum of the absorbed light lost
        backi;%sum of the back scattered light lost
        cutoffi;%sum of the untracked light loss
        housi;%sum of the light lost by absorbing into the housing/attenuating in air before returning 
        approxi;%sum of light that could have been followed, but was assumed to not factor in anymore
        measi;%sum of the measured light
        transi;%sum of transmitted light
        b2hi;%sum of light returning to housing
    end
end