function SystemParam = setupSystemParams()
    SystemParam = SysParam;

SystemParam.division=1*10^4;%2 cm division of the measurements of the figure
%light variables
SystemParam.uvWavelength = 280;    % wavelength of UV light in nmisplacement dependent radiation for LED model
SystemParam.c= 299792458;%[m/s] speed of light
SystemParam.angularFreq=2*pi*SystemParam.c/(SystemParam.uvWavelength*10^-9);%[rad/s]


%parameters for properties of medium a light ray is traveling through
SystemParam.waterInterface = 1;     %is the fiber submerged? no? then it's in air
SystemParam.waterStart=4.5*10^4;      %location at which the water starts (if water is there) 4.5cm [um]
SystemParam.isFiberLossy=1;%if the fiber is considered lossy
if SystemParam.isFiberLossy==1
    SystemParam.n1 = 1.5+1.5E-07*1i;%1.49+3.8E-07;%for PMMA;%%(0.00002317*1i);%1*10^-5*1i);              % RI of Quartz Optical Fiber
else%if not considered a lossy medium then it has no imaginary component
    SystemParam.n1 = 1.50;
end
SystemParam.isMediumLossy=1;%if the medium is considered lossy
SystemParam.kAir = 3*10^-5;             % Attenuation Constant of air 1/cm (assuming 20degC) https://www.ndt.net/article/ultragarsas/63-2008-no.1_03-jakevicius.pdf
SystemParam.kWater =0.0293; % 0.24;% zhe feedwater absorbance at 275nm/cm  %0.0293;         %Attenuation Constant of water 1/cm: https://www.sciencedirect.com/science/article/pii/1350448795002847
SystemParam.kCytop=2.5;%attenuation constant of cytop [cm-1], technical info from cytop page says 5% absorbance over 200um -> 0.05/0.02cm=2.5 cm-1
%calculate the imaginary parts of the RI for air and water from the
%refractive indexes and attenuation coefficients
if SystemParam.isMediumLossy==1
    i_nair=(SystemParam.kAir*10^2)*SystemParam.c/(2*SystemParam.angularFreq);
    i_nWater=(SystemParam.kWater*10^2)*SystemParam.c/(2*SystemParam.angularFreq);

    SystemParam.n2 = 1.00029477+1i*i_nair;              %(Ciddor, 1996) RI of medium w/in the coupling distance/separation distance/"gap" of LED and fiber face
    SystemParam.nWater = 1.353+1i*i_nWater;         % (Hale & Querry,1973) RI of water
else%if not considered a lossy medium then it has no imaginary component
    SystemParam.nWater = 1.353;         % (Hale & Querry,1973) RI of water
    SystemParam.n2 = 1.00029477;              % (Ciddor, 1996) RI of medium w/in the coupling distance/separation distance/"gap" of LED and fiber face
end
SystemParam.nMetal=(1.39+2.3i);%(2.3*1i);%1+(1i*2);% from https://pubs-aip-org.ezproxy1.lib.asu.edu/aip/jap/article/53/9/6340/308961/Optical-constants-and-spectral-selectivity-of

if (SystemParam.waterInterface == 0)
    SystemParam.n5 = SystemParam.n2;          % RI of air
    SystemParam.k = SystemParam.kAir;   % Attenuation Constant of air at 250nm cm^-1 (assuming exceptionally clear air )https://thesis.library.caltech.edu/3249/1/Baum_wa_1950.pdf
else
    SystemParam.n5 = SystemParam.nWater;         % RI of water
    SystemParam.k = SystemParam.kWater;         %Attenuation Constant of water 1/cm: https://www.sciencedirect.com/science/article/pii/1350448795002847
end
%paramaters of the LED
SystemParam.initialIntensity = 100*10^3;      % LED intensity, 100 mW->uW
SystemParam.ledDiameter=1.1;                 %mm, diameter of LED
SystemParam.ledAngle=deg2rad(175);      %angle distribution of light
SystemParam.q=0;                    %exponent for vertical DOM for LED from /science/article/pii/1350448795002847
SystemParam.ledHalfAngle=deg2rad(70);     %Angle from LED manufacturer's at which the power is <1/2 the max
SystemParam.ledDistance = 1*10^3;       %distance of LED from fiber in mm->um
SystemParam.SMA=1;%is there a SMA connector used to couple the fiber to the LED? 1 for yes, 0 for no
SystemParam.smaFlushLength=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SystemParam.smaTotalLength=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SystemParam.smaDiameter=2.5*10^3;%diameter of SMA connector [mm to um]
SystemParam.smaFillLength=0;%2.5*10^4;      %depth the cytop fills the sma connector
SystemParam.isSmaSealed=0;%1;                       %is the SMA filled with cytop?
SystemParam.metalAbsorbPct=0.0;         %absorption fom metal of the SMA connector
%parameters of the optical fiber
SystemParam.fiberRadius = 500/2;       %um
SystemParam.xLen = (12.5*10^4);      %10cm given in (um), length along fiber
SystemParam.measDistance = 0;    %distance measured from the fiber
SystemParam.numFibers=1;            %number of fibers, 1, 4 , or 19
%parameters for changing fidelity of the simulation
%TODO: WAS 21
SystemParam.angleNum=4;               %number of angles of led emission (must be less than numLedRays, below)
SystemParam.numLedRays=SystemParam.angleNum^2;  %total number of rays of led emission
%TODO: WAS 15
SystemParam.scatterNum=8;             %maximum number of rays produced by the scattering during end reflect
SystemParam.frontScatter=1;                       %make this a system param if beneficial to keep it %if there's a scatter cone considered in the front
SystemParam.maxBounce=1;%maximum number of bounces to track
SystemParam.contDx=10000;%continuous transmission interval [um]
SystemParam.maxScatterAngle=pi/18;%3*pi/7;%3*pi/7;%maximum scattering angle 5 degrees
SystemParam.housingBounce=1;%maximum number of times I'm willing to let the ray bounce between the SMA and the fiber
SystemParam.angleDiv=2;%divisffions of angles to look at
SystemParam.difTolerance=10^-10;%tolerance of the difference
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
SystemParam.minPhotons=6.626*(10^-34)*(3*10^8/(SystemParam.uvWavelength*10^-9));%minimum  photon energy in Joules/photon, assumign 1 photon/s
SystemParam.intensityMin=SystemParam.minPhotons*10^15;%(uW)

%absorption/attenuation coefficient calculations, Gerd Keiser. (2011). Optical Fiber Communications (4th ed.). McGraw-Hill Education.
%attenuation coefficients
SystemParam.rayScatterCoeff=0.9;%0.125;%,0.9]; %coefficient for what percentage of the loss is attributed to scattering vs absorption

% 
% friom https://www.content.molex.com/dxdam/literature/987650-8936.pdf
%SystemParam.alpha_Cytop=(1000*10^-9);%db/um from cytop estimate in friday pres 091622
%38.9 38.9 3968*10^-6;%197.55,98.33*10^-6;;%12.99*10^-6;%db/um%(1.64*10^-9)*(850/SystemParam.uvWavelength)^4;%db/um 
%SystemParam.rayleighCoeff=(1-10^(-SystemParam.rayleighAlpha*SystemParam.xLen/10));%coefficient used to find the amt of light scattered
%use function to calculate the percent side scattered in the theta_OC
[Scat_distrib] = Scatter_Coeff(SystemParam);
SystemParam.scatterCoeff=Scat_distrib;