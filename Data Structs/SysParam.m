classdef SysParam
    
    properties
        n1;                 % RI of QOF
        n2;                 % RI of air (inside the coupling gap)
        n3;                 % RI of NP
        n4;                 % RI of polymer cytop
        n5;                 % RI of external medium, either air or water
        nwater;             %specific RI for water
        xlen;               % length of fiber, in um
        meas_distance;      % distance away from fiber, in um
       
        k_water;            % UV atten coeff of water
        kair;               % attenuation constant of air
        k;                  % general attenuation constant of media
        r_fiber;            % radius of fiber, in um
        uv_wavelength;      % wavelength of UV light, in nm
        
        led_distance;       % distance between LED and fiber
        I_init;             % intensity of LED
        poly_coat;          % boolean if a polymer coating was applied
        np_coat;            % boolean if a polymer coating was applied
        LEDd;               % diameter of LED
        LEDang;             % full angle of radiation of LED
        LEDa;               % full angle of radiation of LED (check why two of these)
        q;                  % parameter for y direction LED intensity changes (possibly unused)
        NPstart;            % where the NP along the fiber starts
        Half_t;             % half angle of the LED
        alpha_r;      %coefficient used in determining amount of rayleigh scattering in glass fiber
        RayleighCoeff;      %coefficient used to determine total loss for entire distance of the fiber
        alpha_uv;           %coefficient used in determining amount of UV absorption occurs in glass fiber
        num_fib;            %numble of fibers in bundle
        photon_min;         %energy per photon in Joules
        raynum;             %number of rays from led emission
        angnum;             %number of angles looked at from led emission
        waterstart;         %location at which water interface starts
        waterInterface;     %if there is a water interface
        metal_abs;          %percentage of light absorbed by SMA connector metal
        division;           %length of division for data display
        scatnum;            %number of scatter rays produced by end reflect
        maxbounce;          %number of times the ray is allowed to bounce from one transmission end to another
        I_min;              %minimum Intensity worthwhile to keep track of
        cont_dx;            %continuous travel dx min
        scat_ang_max;       %maximum scattering angle
        dxordr;             %logical selection whether using dx(=0) or dr(=1) for attenuation changes
        housing_bounce;     %maximum number of bounces considered in the sma housing
        ang_div;            %how many angles will be considered
        SMA;
        SMA_flushlength;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
        SMA_totallength;%the total length of the SMA connector is 2.5 cm long
        SMA_diam;
        scatter_coeff;
        dif_tol;
        
        %added on 8/31/23
        lossy_f=1;%if the fiber is considered lossy
        lossy_m=1;%if the medium is considered lossy
        n_metal;%refractive index of the metal as estimated from https://pubs-aip-org.ezproxy1.lib.asu.edu/aip/jap/article/53/9/6340/308961/Optical-constants-and-spectral-selectivity-of
        c;%speed of light
        ang_freq;%angular frequency
    end
end
