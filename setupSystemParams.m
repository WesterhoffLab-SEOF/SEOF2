function SystemParam = setupSystemParams()
    SystemParam = SysParam;

SystemParam.division=1*10^4;%2 cm division of the measurements of the figure
%light variables
SystemParam.uv_wavelength = 280;    % wavelength of UV light in nmisplacement dependent radiation for LED model
SystemParam.c= 299792458;%[m/s] speed of light
SystemParam.ang_freq=2*pi*SystemParam.c/(SystemParam.uv_wavelength*10^-9);%[rad/s]


%parameters for properties of medium a light ray is traveling through
SystemParam.waterInterface = 1;     %is the fiber submerged? no? then it's in air
SystemParam.waterstart=4.5*10^4;      %location at which the water starts (if water is there) 4.5cm [um]
SystemParam.lossy_f=1;%if the fiber is considered lossy
if SystemParam.lossy_f==1
    SystemParam.n1 = 1.5+1.5E-07*1i;%1.49+3.8E-07;%for PMMA;%%(0.00002317*1i);%1*10^-5*1i);              % RI of Quartz Optical Fiber
else%if not considered a lossy medium then it has no imaginary component
    SystemParam.n1 = 1.50;
end
SystemParam.lossy_m=1;%if the medium is considered lossy
SystemParam.kair = 3*10^-5;             % Attenuation Constant of air 1/cm (assuming 20degC) https://www.ndt.net/article/ultragarsas/63-2008-no.1_03-jakevicius.pdf
SystemParam.k_water =0.0293; % 0.24;% zhe feedwater absorbance at 275nm/cm  %0.0293;         %Attenuation Constant of water 1/cm: https://www.sciencedirect.com/science/article/pii/1350448795002847
SystemParam.k_cytop=2.5;%attenuation constant of cytop [cm-1], technical info from cytop page says 5% absorbance over 200um -> 0.05/0.02cm=2.5 cm-1
%calculate the imaginary parts of the RI for air and water from the
%refractive indexes and attenuation coefficients
if SystemParam.lossy_m==1
    i_nair=(SystemParam.kair*10^2)*SystemParam.c/(2*SystemParam.ang_freq);
    i_nwater=(SystemParam.k_water*10^2)*SystemParam.c/(2*SystemParam.ang_freq);

    SystemParam.n2 = 1.00029477+1i*i_nair;              %(Ciddor, 1996) RI of medium w/in the coupling distance/separation distance/"gap" of LED and fiber face
    SystemParam.nwater = 1.353+1i*i_nwater;         % (Hale & Querry,1973) RI of water
else%if not considered a lossy medium then it has no imaginary component
    SystemParam.nwater = 1.353;         % (Hale & Querry,1973) RI of water
    SystemParam.n2 = 1.00029477;              % (Ciddor, 1996) RI of medium w/in the coupling distance/separation distance/"gap" of LED and fiber face
end
SystemParam.n_metal=(1.39+2.3i);%(2.3*1i);%1+(1i*2);% from https://pubs-aip-org.ezproxy1.lib.asu.edu/aip/jap/article/53/9/6340/308961/Optical-constants-and-spectral-selectivity-of

if (SystemParam.waterInterface == 0)
    SystemParam.n5 = SystemParam.n2;          % RI of air
    SystemParam.k = SystemParam.kair;   % Attenuation Constant of air at 250nm cm^-1 (assuming exceptionally clear air )https://thesis.library.caltech.edu/3249/1/Baum_wa_1950.pdf
else
    SystemParam.n5 = SystemParam.nwater;         % RI of water
    SystemParam.k = SystemParam.k_water;         %Attenuation Constant of water 1/cm: https://www.sciencedirect.com/science/article/pii/1350448795002847
end
%paramaters of the LED
SystemParam.I_init = 100*10^3;      % LED intensity, 100 mW->uW
SystemParam.LEDd=1.1;                 %mm, diameter of LED
SystemParam.LEDa=deg2rad(175);      %angle distribution of light
SystemParam.q=0;                    %exponent for vertical DOM for LED from /science/article/pii/1350448795002847
SystemParam.Half_t=deg2rad(70);     %Angle from LED manufacturer's at which the power is <1/2 the max
SystemParam.led_distance = 1*10^3;       %distance of LED from fiber in mm->um
SystemParam.SMA=1;%is there a SMA connector used to couple the fiber to the LED? 1 for yes, 0 for no
SystemParam.SMA_flushlength=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SystemParam.SMA_totallength=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SystemParam.SMA_diam=2.5*10^3;%diameter of SMA connector [mm to um]
SystemParam.SMA_fill=0;%2.5*10^4;      %depth the cytop fills the sma connector
SystemParam.sealed_SMA=0;%1;                       %is the SMA filled with cytop?
SystemParam.metal_abs=0.0;         %absorption fom metal of the SMA connector
%parameters of the optical fiber
SystemParam.r_fiber = 500/2;       %um
SystemParam.xlen = (12.5*10^4);      %10cm given in (um), length along fiber
SystemParam.meas_distance = 0;    %distance measured from the fiber
SystemParam.num_fib=1;            %number of fibers, 1, 4 , or 19
%parameters for changing fidelity of the simulation
%TODO: WAS 21
SystemParam.angnum=4;               %number of angles of led emission (must be less than raynum, below)
SystemParam.raynum=SystemParam.angnum^2;  %total number of rays of led emission
%TODO: WAS 15
SystemParam.scatnum=8;             %maximum number of rays produced by the scattering during end reflect
SystemParam.scat_front=1;                       %make this a system param if beneficial to keep it %if there's a scatter cone considered in the front
SystemParam.maxbounce=1;%maximum number of bounces to track
SystemParam.cont_dx=10000;%continuous transmission interval [um]
SystemParam.scat_ang_max=pi/18;%3*pi/7;%3*pi/7;%maximum scattering angle 5 degrees
SystemParam.housing_bounce=1;%maximum number of times I'm willing to let the ray bounce between the SMA and the fiber
SystemParam.ang_div=2;%divisffions of angles to look at
SystemParam.dif_tol=10^-10;%tolerance of the difference
%possibly include SystemParam for point at which a ray is considered
%horizontal or vertical
%trying out different logic for stuff
%SystemParam.dxordr=1;%select 0 if the logic is using dX, select 1 if the logic is using the whole travel distance

%include parameters for the polymer coating factors
SystemParam.n4=1.353;
%include parameters for the nanoparticle coating factors
SystemParam.n3=1.5; %RI of silica beads

%calculated parameters
%minimum photon energy at the wavelength of light before it changes "color"
SystemParam.photon_min=6.626*(10^-34)*(3*10^8/(SystemParam.uv_wavelength*10^-9));%minimum  photon energy in Joules/photon, assumign 1 photon/s
SystemParam.I_min=SystemParam.photon_min*10^15;%(uW)

%absorption/attenuation coefficient calculations, Gerd Keiser. (2011). Optical Fiber Communications (4th ed.). McGraw-Hill Education.
%attenuation coefficients
SystemParam.r_coeff=0.9;%0.125;%,0.9]; %coefficient for what percentage of the loss is attributed to scattering vs absorption

% 
% friom https://www.content.molex.com/dxdam/literature/987650-8936.pdf
%SystemParam.alpha_Cytop=(1000*10^-9);%db/um from cytop estimate in friday pres 091622
%38.9 38.9 3968*10^-6;%197.55,98.33*10^-6;;%12.99*10^-6;%db/um%(1.64*10^-9)*(850/SystemParam.uv_wavelength)^4;%db/um 
%SystemParam.RayleighCoeff=(1-10^(-SystemParam.alpha_r*SystemParam.xlen/10));%coefficient used to find the amt of light scattered
%use function to calculate the percent side scattered in the theta_OC
[Scat_distrib] = Scatter_Coeff(SystemParam);
SystemParam.scatter_coeff=Scat_distrib;