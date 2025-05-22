classdef SysParam
    
    properties
        n1;                 % RI of QOF
        n2;                 % RI of air (inside the coupling gap)
        n3;                 % RI of NP
        n4;                 % RI of polymer cytop
        n5;                 % RI of external medium, either air or water
        nWater;             %specific RI for water
        xLen;               % length of fiber, in um
        measDistance;      % distance away from fiber, in um
       
        kWater;            % UV atten coeff of water
        kAir;               % attenuation constant of air
        k;                  % general attenuation constant of media
        fiberRadius;            % radius of fiber, in um
        uvWavelength;      % wavelength of UV light, in nm
        
        ledDistance;       % distance between LED and fiber
        initialIntensity;             % intensity of LED
        polymerCoat;          % boolean if a polymer coating was applied
        nanoParticleCoat;            % boolean if a polymer coating was applied
        ledDiameter;               % diameter of LED
        ledAngle;             % full angle of radiation of LED
        q;                  % parameter for y direction LED intensity changes (possibly unused)
        npStartingPoint;            % where the NP along the fiber starts
        ledHalfAngle;             % half angle of the LED
        rayleighAlpha;      %coefficient used in determining amount of rayleigh scattering in glass fiber
        rayleighCoeff;      %coefficient used to determine total loss for entire distance of the fiber
        uvAlpha;           %coefficient used in determining amount of UV absorption occurs in glass fiber
        alpha;              %total loss coefficient of the optical fiber
        rayScatterCoeff;            %proportion of total loss going to rayleigh scattering
        numFibers;            %numble of fibers in bundle
        minPhotons;         %energy per photon in Joules
        numLedRays;             %number of rays from led emission
        angleNum;             %number of angles looked at from led emission
        waterStart;         %location at which water interface starts
        waterInterface;     %if there is a water interface
        isSmaSealed;         %is the SMA filled w cytop
        metalAbsorbPct;          %percentage of light absorbed by SMA connector metal
        division;           %length of division for data display
        scatterNum;            %number of scatter rays produced by end reflect
        maxBounce;          %number of times the ray is allowed to bounce from one transmission end to another
        intensityMin;              %minimum Intensity worthwhile to keep track of
        contDx;            %continuous travel dx min
        maxScatterAngle;       %maximum scattering angle
        dxordr;             %logical selection whether using dx(=0) or dr(=1) for attenuation changes
        housingBounce;     %maximum number of bounces considered in the sma housing
        angleDiv;            %how many angles will be considered
        SMA;
        smaFlushLength;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
        smaTotalLength;%the total length of the SMA connector is 2.5 cm long
        smaDiameter;
        scatterCoeff;
        difTolerance;
        frontScatter;
        
        %added on 8/31/23
        isFiberLossy=1;%if the fiber is considered lossy
        isMediumLossy=1;%if the medium is considered lossy
        nMetal;%refractive index of the metal as estimated from https://pubs-aip-org.ezproxy1.lib.asu.edu/aip/jap/article/53/9/6340/308961/Optical-constants-and-spectral-selectivity-of
        c;%speed of light
        angularFreq;%angular frequency
        kCytop;%attenuation of the cytop
        smaFillLength; %length the sma connector is filled to
    end
end
