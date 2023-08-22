%%simplest possible fiber model%%
%for the most part, the comments that use the following phrases mean:
%to do =  i need to do it before the model will run properly
%change =  I should wait to change a variable name until I can do it
%en-masse
%include =  its something for later when the model has grown
close all; clear all;%fresh slate
addpath('Utilities');%functions go into utilities folder
addpath('Data Structs');%structs go into this folder

%%%%%%%%%%%%SETTING UP THE SIMULATION%%%%%%%%%%%%%%
%systemParam, direct input, can change these
SystemParam = SysParam;
%data display variables... what do you want to see generated?
writeToFile =1;   % set to 1 to save generated results to excel file, set to 0 otherwise
filename = 'C:\Users\Nora Shapiro\OneDrive - Arizona State University\SEOF_RayTracing_071823\data\MODEL_simple_081723.xlsx'; %name of file %MUST INCLUDE ADDRESS OF DESTINATION FOLDER    
sheet_num=4;%Excel Sheet number 
description="Changes: trying a new ray model out. doing 5x5 ray, 5 scatter, 1 bounce max, 2bounce sma . yes sma";%description of the changes made in this file
SystemParam.division=1*10^4;%2 cm division of the measurements of the figure

%parameters for properties of medium a light ray is traveling through
SystemParam.waterInterface = 0;     %is the fiber submerged? no? then it's in air
SystemParam.waterstart=1*10^4;      %location at which the water starts (if water is there) 4.5cm [um]
SystemParam.n1 = 1.50;              % RI of Quartz Optical Fiber
SystemParam.n2 = 1.00;              % RI of medium w/in the coupling distance/separation distance/"gap" of LED and fiber face
SystemParam.nwater = 1.333;         % RI of water
SystemParam.kair = 3*10^-5;             % Attenuation Constant of air 1/cm (assuming 20degC) https://www.ndt.net/article/ultragarsas/63-2008-no.1_03-jakevicius.pdf
if (SystemParam.waterInterface == 0)
    ext_media="External Medium: Air"; %string for excel document description
    SystemParam.n5 = 1.00;          % RI of air
    SystemParam.k = 3/(1*10^5);   % Attenuation Constant of air at 250nm cm^-1 (assuming exceptionally clear air )https://thesis.library.caltech.edu/3249/1/Baum_wa_1950.pdf
else
    SystemParam.n5 = 1.333;         % RI of water
    SystemParam.k = 0.0293;         %Attenuation Constant of water 1/cm: https://www.sciencedirect.com/science/article/pii/1350448795002847
    ext_media="External Medium: Air"; %string for excel document description
end
%paramaters of the LED
SystemParam.I_init = 70*(10^3);      % LED intensity, 70 mW
SystemParam.LEDd=0.8;                 %mm, diameter of LED
SystemParam.LEDa=deg2rad(175);      %angle distribution of light
SystemParam.q=0;                   %exponent for vertical DOM for LED from /science/article/pii/1350448795002847
SystemParam.Half_t=deg2rad(70);     %Angle from LED manufacturer's at which the power is <1/2 the max
SystemParam.led_distance = 1;       %distance of LED from fiber in mm
SystemParam.uv_wavelength = 275;    % wavelength of UV light in nmisplacement dependent radiation for LED model
SystemParam.SMA=1;%is there a SMA connector used to couple the fiber to the LED? 1 for yes, 0 for no
SystemParam.SMA_flushlength=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SystemParam.SMA_totallength=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SystemParam.SMA_diam=2.5*10^3;%diameter of SMA connector [mm to um]
SystemParam.metal_abs=0.0;         %absorption from metal of the SMA connector
%parameters of the optical fiber
SystemParam.r_fiber = 500/2;       %um
SystemParam.xlen = (12.5*10^4);      %10cm given in (um), length along fiber
SystemParam.meas_distance = 0;    %distance measured from the fiber
SystemParam.num_fib=1;                %number of fibers, 1, 4 , or 19
%parameters for changing fidelity of the simulation
SystemParam.angnum=5;               %number of angles of led emission (must be less than raynum, below)
SystemParam.raynum=SystemParam.angnum^2;  %total number of rays of led emission
SystemParam.scatnum=5;%maximum number of rays produced by the scattering during end reflect
SystemParam.maxbounce=1;%100;%maximum number of bounces to track
SystemParam.cont_dx=2*10^4;%continuous transmission interval [um]
SystemParam.scat_ang_max=pi/3;%maximum scattering angle
SystemParam.housing_bounce=2;%maximum number of times I'm willing to let the ray bounce between the SMA and the fiber
SystemParam.ang_div=2;%divisions of angles to look at
SystemParam.dif_tol=10^-6;%tolerance of the difference
%possibly include SystemParam for point at which a ray is considered
%horizontal or vertical
%trying out different logic for stuff
%SystemParam.dxordr=1;%select 0 if the logic is using dX, select 1 if the logic is using the whole travel distance

%include parameters for the polymer coating factors

%include parameters for the nanoparticle coating factors

%calculated parameters
%minimum photon energy at the wavelength of light before it changes "color"
SystemParam.photon_min=6.626*(10^-34)*(3*10^8/(SystemParam.uv_wavelength*10^-9));%minimum  photon energy in Joules/photon, assumign 1 photon/s
SystemParam.I_min=SystemParam.photon_min*10;%(uW)

%absorption/attenuation coefficient calculations, Gerd Keiser. (2011). Optical Fiber Communications (4th ed.). McGraw-Hill Education.
%attenuation coefficients
SystemParam.alpha_uv=(0.1*10^-6);%db/um friom https://www.content.molex.com/dxdam/literature/987650-8936.pdf
%SystemParam.alpha_Cytop=(1000*10^-9);%db/um from cytop estimate in friday pres 091622
SystemParam.alpha_r=(1.64*10^-9)*(850/SystemParam.uv_wavelength)^4;%db/um 
%SystemParam.RayleighCoeff=(1-10^(-SystemParam.alpha_r*SystemParam.xlen/10));%coefficient used to find the amt of light scattered
%use function to calculate the percent side scattered in the theta_OC
[Scat_distrib] = Scatter_Coeff(SystemParam);
SystemParam.scatter_coeff=Scat_distrib;

%%%SETTING UP PARAMETERS NEEDED FOR ITERATION AND DATA DISPLAY%%%%
%all of these variables can be turned into a vector of different values so
%it will automatically re-run the model with a new variable.
%CAN ONLY ITERATE ON ONE AT A TIME FOR NOW
r_fiber=SystemParam.r_fiber;%[125,250,500]%um alternate
xlen=SystemParam.xlen;%([10,15,20,25,30,35,40,45,50]+2)*10^4;%cm->[um] alternate
led_dist=SystemParam.led_distance;% [0,1,2,3,4,5];%mm alternate
%include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vector of the size of each possible iterateable vector
%lengths of each iterable variable
lengths_it=[length(r_fiber),length(xlen),length(led_dist)];%include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it_num=max(lengths_it);
index_it=find(lengths_it>1);%indexes of iterable variable(s)

%if/else statement for use in displaying data about the iteration/fiber
%bundle
if ~isempty(index_it)%if there are iterable variables
    if length(index_it)>1%if there are multiple variables with multiple iterations, end and set up again
        yes_itname=0;
        error('More than one variable is being iterated, please fix and re-try')
    elseif SystemParam.num_fib>1 %if there is more than one fiber in the bundle and we're iterating (this should be a capability we have eventually though)
        yes_itname=0;
        error('You are trying to iterate different variables for a bundle of fibers, please fix and re-try')
    else%only one variable with multiple iteration. Determine which variable is being iterated
        yes_itname=1;%name of the iteration needs to be incorporated into data display
        switch index_it
            case 1
                it_name='Radius of Fiber [/mum]';
                Title_Main=sprintf('Irradiance along the Fiber Length, Varying Fiber Diameter \n %d Fiber(s), %d cm long, %d mm coupling distance',SystemParam.num_fib, xlen/10^4,led_dist);
                leg_Main=cell{lengths_it(index_it)};
                    for i=1:lengths_it(index_it)
                        leg_Main(i)={sprintf('%d [\mum] Diameter',r_fiber(i)*2)};
                    end
                legend_Main=string(leg_Main);
            case 2
                it_name='Length of Fiber [/mum]';
                Title_main=sprintf('Irradiance along the Fiber Length, Varying Fiber Length \n %d Fiber(s),  %d [um] OD , %d mm coupling distance',SystemParam.num_fib,r_fiber*2,led_dist);
                leg_Main=cell{lengths_it(index_it)};
                    for i=1:lengths_it(index_it)
                        leg_Main(i)={sprintf('%d [cm] Length',xlen(i)/10^4)};
                    end
                legend_Main=string(leg_Main);
            case 3
                it_name='LED Distance from Fiber [mm]';
                Title_Main=sprintf('Irradiance along the Fiber Length, Varying Coupling Distance \n %d Fiber(s), %d cm long, %d [um] OD',SystemParam.num_fib, xlen/10^4,r_fiber*2);
                leg_Main=cell{lengths_it(index_it)};
                    for i=1:lengths_it(index_it)
                        leg_Main(i)={sprintf('%d [mm] Coupling Distance',led_dist(i))};
                    end
                legend_Main=string(leg_Main);
                %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            otherwise 
                error('unexpected iteratable variable. Please add a new switch-case for it')
        end
    end

else%if there are no iterable variables (just running 1 case)
    if SystemParam.num_fib>1%if there are no iterable variables but there are a bundle of fibers to look at
    it_name='Fiber in Bundle';
    yes_itname=1;%iteration name needed for data display
    leg_Main=cell(1,SystemParam.num_fib);%initialize cell for string vector for the legend
    for i=1:SystemParam.num_fib
        leg_Main(i)={sprintf('Fiber %d',i)};
    end
    legend_Main=string(leg_Main);
    else %only one fiber, don't need a legend
        yes_itname=0;%no iteration name needed for data display
        it_name=' ';
    end
    Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,r_fiber*2,led_dist);
end

%set up number of random smoothing loops
 aa_lim=1;%include if/then statements for increasing the limit to higher if there are random input parameters
 %Maybe there should be a system param for how smoothed the data should be?

%num_max = floor((SystemParam.xlen/10e3)/(SysteParam.division/10e3));%dividing length of fiber by the increments, then adding one to have data at each end of bin
%accessing data storage structs
ray=Ray; %for each individual ray
aa_Rand=aaRand; %for each trial run to smooth out random values
FibIt=Fib_It;%for each fiber
SumVal=SummaryVals;%ALL of the ray data for each fiber. cell.
Tally=Tracking;%tracking how much certain things happen. Assists with troubleshooting
meas=Measure_Instant;%instant, temporary storage of measurable data
IT=travel_storage;%temporary storage of light losses within the travel function
%empty vectors for storing all of the summary data of each fiber at each
%iteration. (1) is the average val, (2) is the std of the avg val
FibIt(1).Pow_enter=zeros(it_num,SystemParam.num_fib);FibIt(2).Pow_enter=zeros(it_num,SystemParam.num_fib);
FibIt(1).transmitted=zeros(it_num,SystemParam.num_fib);FibIt(2).transmitted=zeros(it_num,SystemParam.num_fib);
FibIt(1).pow_side=zeros(it_num,SystemParam.num_fib);FibIt(2).pow_side=zeros(it_num,SystemParam.num_fib);
FibIt(1).pow_side_use=zeros(it_num,SystemParam.num_fib);FibIt(2).pow_side_use=zeros(it_num,SystemParam.num_fib);
FibIt(1).SMA_pow_side=zeros(it_num,SystemParam.num_fib);FibIt(2).SMA_pow_side=zeros(it_num,SystemParam.num_fib);
FibIt(1).pow_side_waterst=zeros(it_num,SystemParam.num_fib);FibIt(2).pow_side_waterst=zeros(it_num,SystemParam.num_fib);
FibIt(1).absorbed=zeros(it_num,SystemParam.num_fib);FibIt(2).absorbed=zeros(it_num,SystemParam.num_fib);
FibIt(1).backscat=zeros(it_num,SystemParam.num_fib);FibIt(2).backscat=zeros(it_num,SystemParam.num_fib);
FibIt(1).SMAabs=zeros(it_num,SystemParam.num_fib);FibIt(2).SMAabs=zeros(it_num,SystemParam.num_fib);
FibIt(1).approxpow_dif=zeros(it_num,SystemParam.num_fib);FibIt(2).approxpow_dif=zeros(it_num,SystemParam.num_fib);
FibIt(1).approxpow_pos=zeros(it_num,SystemParam.num_fib);FibIt(2).approxpow_pos=zeros(it_num,SystemParam.num_fib);
FibIt(1).b2hpow=zeros(it_num,SystemParam.num_fib);FibIt(2).b2hpow=zeros(it_num,SystemParam.num_fib);
FibIt(1).cutoffpow=zeros(it_num,SystemParam.num_fib);FibIt(2).cutoffpow=zeros(it_num,SystemParam.num_fib);
FibIt(1).remaininglosses=zeros(it_num,SystemParam.num_fib);FibIt(2).remaininglosses=zeros(it_num,SystemParam.num_fib);
FibIt(1).UC=zeros(it_num,SystemParam.num_fib);FibIt(2).UC=zeros(it_num,SystemParam.num_fib);
FibIt(1).RatioIsIt=zeros(it_num,SystemParam.num_fib);FibIt(2).RatioIsIt=zeros(it_num,SystemParam.num_fib);
FibIt(1).Y=cell(it_num,SystemParam.num_fib);FibIt(2).Y=cell(it_num,SystemParam.num_fib);

%empty cell setup for each individual ray tracking info
SumVal.Pow_enter=cell(it_num,SystemParam.num_fib);
SumVal.transmitted=cell(it_num,SystemParam.num_fib);
SumVal.pow_side=cell(it_num,SystemParam.num_fib);
SumVal.pow_side_use=cell(it_num,SystemParam.num_fib);
SumVal.SMA_pow_side=cell(it_num,SystemParam.num_fib);
SumVal.pow_side_waterst=cell(it_num,SystemParam.num_fib);
SumVal.absorbed=cell(it_num,SystemParam.num_fib);
SumVal.backscat=cell(it_num,SystemParam.num_fib);
SumVal.SMAabs=cell(it_num,SystemParam.num_fib);
SumVal.approxpow_dif=cell(it_num,SystemParam.num_fib);%sum of power differences due to difference in approximations
SumVal.approxpow_pos=cell(it_num,SystemParam.num_fib); %absolute value of all of the power difference approximationsSumVal.b2hpow=cell(it_num,SystemParam.num_fib);
SumVal.cutoffpow=cell(it_num,SystemParam.num_fib);
SumVal.remaininglosses=cell(it_num,SystemParam.num_fib);
SumVal.UC=cell(it_num,SystemParam.num_fib);
SumVal.RatioIsIt=cell(it_num,SystemParam.num_fib);
SumVal.Y=cell(it_num,SystemParam.num_fib);

%%%%%%EACH ITERATION OF THE FIBER
for iteration=1:it_num
 
    %%%%%%%reset parameters to be reflective of the iteration number
    if yes_itname==1%if there are iterable items
        if index_it==1
            SystemParam.r_fiber=r_fiber(iteration);
        elseif index_it==2
            SystemParam.xlen=xlen(iteration);
        elseif index_it==3
            SystemParam.led_distance=led_dist(iteration);
        %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
Xvec=0:(SystemParam.division/10e3):(SystemParam.xlen/10e3);
%%%%%%%setting up information for req number of random loops%%%%
    %include if/then surface roughness/NP randomness and such statements
    %statements%%%%%%%%%%%%
    
    %%%%% create the boundary struct for a glass fiber%%%%%%%
    Bounds = Bound;
    Bounds.Pf0=[0,SystemParam.r_fiber];%boundary at starting position
    Bounds.Pfe=[SystemParam.xlen,SystemParam.r_fiber];%boundary at ending position
    
    %%Uniformity coefficient index set up for later calculation
    i1=round((((xlen/10000)-1)/(SystemParam.division*10^-4))*(0.1))+1;
    i6=round((((xlen/10000)-2)/(SystemParam.division*10^-4))*(0.6))+1;
    
    
    %set up of incoming ray angle, intensity, for a variable number of fibers
    LED_d=SystemParam.led_distance*(10^-3);%(mm)%LED distance in mm for Entering Light
    [Ray_X,Ray_Y,alpha_ang,beta_ang,Power_mat] = OutputLED3D(SystemParam);
    powercheck=sum(Power_mat,'all');% make sure power isn't greater than the total LED power
    %determining light rays entering each fiber
    [circles,combocircles] = fiberBundle3D(SystemParam.num_fib,r_fiber*2);
    %[Ent_Int,Theta,y0,perc_hit,perc_ent,Entering_Int,Entering_angle,Entering_X,Entering_Y,incoming_int,incoming_ang] = enterLight3D_AINsubstrate(SystemParam,r_fiber, LED_d,alpha_ang,beta_ang,Ray_X,Ray_Y,Intensity_mat,circles);
    [Ent_Pow,Theta,y0,perc_hit,perc_ent,Entering_Int,Entering_angle,Entering_X,Entering_Y,incoming_int,incoming_ang] = enterLight3D(SystemParam,r_fiber, LED_d,alpha_ang,beta_ang,Ray_X,Ray_Y,Power_mat,circles);

    [a,b,c]=size(Ent_Pow); %[num ang changes,num y locations, number fibers]
              %empty storage for a numeric count of some stuff for trouble
              %shooting
        Tally.photon_vert_count=zeros(it_num,SystemParam.num_fib,aa_lim);
        Tally.photon_min_count=zeros(it_num,SystemParam.num_fib,aa_lim);
        Tally.transmission_count=zeros(it_num,SystemParam.num_fib,aa_lim);%all the time the transmission count is reached
        Tally.big_dif_count_main=zeros(it_num,SystemParam.num_fib,aa_lim);%counting all of the times theres a significant difference with a magnitude >10^-6
        Tally.big_dif_amt_main=zeros(it_num,SystemParam.num_fib,aa_lim);%recording the significant difference amount
        Tally.big_dif_pos_main=zeros(it_num,SystemParam.num_fib,aa_lim);%recording the significant difference amount
        Tally.big_dif_fun_main=cell(it_num,SystemParam.num_fib,aa_lim);%record the name of the function causing issues
        Tally.big_dif_count_meas_main=cell(it_num,SystemParam.num_fib,aa_lim);
        Tally.IO_count_main=zeros(it_num,SystemParam.num_fib,aa_lim);%counting all of the times theres a significant difference with a magnitude >10^-6
        Tally.IO_amt_main=zeros(it_num,SystemParam.num_fib,aa_lim);%recording the significant difference amount
        Tally.IO_pos_main=zeros(it_num,SystemParam.num_fib,aa_lim);%recording the significant difference amount
        Tally.IO_fun_main=cell(it_num,SystemParam.num_fib,aa_lim);%record the name of the function causing issues
        Tally.IO_count_meas_main=cell(it_num,SystemParam.num_fib,aa_lim);

        %within a function
        Tally.big_dif_count_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%counting all of the times theres a significant difference with a magnitude >10^-6
        Tally.big_dif_amt_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%recording the significant difference amount
        Tally.big_dif_pos_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%recording the significant difference amount
        Tally.big_dif_fun_func=cell(it_num,SystemParam.num_fib,aa_lim,a,b);%record the name of the function causing issues
        Tally.big_dif_count_meas_func=cell(it_num,SystemParam.num_fib,aa_lim,a,b);
        Tally.IO_count_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%counting all of the times theres a significant difference with a magnitude >10^-6
        Tally.IO_amt_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%recording the significant difference amount
        Tally.IO_pos_func=zeros(it_num,SystemParam.num_fib,aa_lim,a,b);%recording the significant difference amount
        Tally.IO_fun_func=cell(it_num,SystemParam.num_fib,aa_lim,a,b);%record the name of the function causing issues
        Tally.IO_count_meas_func=cell(it_num,SystemParam.num_fib,aa_lim,a,b);
        Tally.case1count=zeros(it_num,SystemParam.num_fib,aa_lim);
        Tally.case2count=zeros(it_num,SystemParam.num_fib,aa_lim);
        Tally.case3count=zeros(it_num,SystemParam.num_fib,aa_lim);
        Tally.whiletravcount=zeros(it_num,SystemParam.num_fib,aa_lim);

    %within each individual fiber
    for h=1:c%%%%%%include for number fiber change
        %empty storage for each aa_lim
        aa_Rand.Pow_enter=zeros(1,aa_lim);
        aa_Rand.transmitted=zeros(1,aa_lim);
        aa_Rand.pow_side=zeros(1,aa_lim);
        aa_Rand.pow_side_use=zeros(1,aa_lim);
        aa_Rand.SMA_pow_side=zeros(1,aa_lim);
        aa_Rand.pow_side_waterst=zeros(1,aa_lim);
        aa_Rand.absorbed=zeros(1,aa_lim);
        aa_Rand.backscat=zeros(1,aa_lim);
        aa_Rand.SMAabs=zeros(1,aa_lim);
        aa_Rand.approxpow=zeros(1,aa_lim);
        aa_Rand.approxpow_dif=zeros(1,aa_lim);
        aa_Rand.approxpow_pos=zeros(1,aa_lim);
        aa_Rand.b2hpow=zeros(1,aa_lim);
        aa_Rand.cutoffpow=zeros(1,aa_lim);
        aa_Rand.remaininglosses=zeros(1,aa_lim);
        aa_Rand(1).UC=zeros(1,aa_lim);aa_Rand(2).UC=zeros(1,aa_lim);
        aa_Rand(1).RatioIsIt=zeros(1,aa_lim);aa_Rand(2).RatioIsIt=zeros(1,aa_lim);
        aa_Rand(1).Y=cell(aa_lim,1);aa_Rand(2).Y=cell(aa_lim,1);
        aa_Rand(1).meas_inten=cell((a*b),aa_lim); %length of all the rays used
        aa_Rand(1).meas_points=cell((a*b),aa_lim);%length of all the rays used
        %empty storage vector setup for each individual ray 
        %possibly will need to include dimension for the scatter cone later
        %(a,b,ggmax,aa_lim)
        ray.Pow_enter=zeros(aa_lim,a,b);%=I_ent0
        ray.transmitted=zeros(aa_lim,a,b);
        ray.pow_side=zeros(aa_lim,a,b);
        ray.pow_side_use=zeros(aa_lim,a,b);
        ray.SMA_pow_side=zeros(aa_lim,a,b);
        ray.pow_side_waterst=zeros(aa_lim,a,b);
        ray.absorbed=zeros(aa_lim,a,b);
        ray.backscat=zeros(aa_lim,a,b);
        ray.SMAabs=zeros(aa_lim,a,b);
        ray.approxpow=zeros(aa_lim,a,b);
        ray.approxpow_dif=zeros(aa_lim,a,b);%sum of power differences due to difference in approximations
        ray.approxpow_pos=zeros(aa_lim,a,b); %absolute value of all of the power difference approximations
        ray.b2hpow=zeros(aa_lim,a,b);
        ray.cutoffpow=zeros(aa_lim,a,b);
         ray.remaininglosses=zeros(aa_lim,a,b);%should theoretically be equal to the approxpow_dif?;
        ray.UC=zeros(aa_lim,a,b);
        ray.RatioIsIt=zeros(aa_lim,a,b);
        ray.Y=cell(aa_lim,a,b);
        ray.meas_inten=cell(aa_lim,a,b);
        ray.meas_points=cell(aa_lim,a,b);
        %tracking and summmation values for each fiber
        for aa=1:aa_lim
            %include surface roughness distribution
            %include NP distribution
                
            for xx=1:a
                for yy=1:b
                %initialize measurement storage vectors
                meas.points= zeros(2*10e6,2);
                meas.inten = zeros(2*10e6,1);
                meas.counter=1;
                meas.sum=0;
                
                Global_Index=[iteration,h,aa,xx,yy];%record the index values, may need to include a gg if scatter cone at the front
                %tracking and summing values for each ray
                %temporary storage of data

                %initialize summed/tallied values
                Y_sum=zeros(size(Xvec));
                
                %include ? scatter cone. would need to go back and add a gg
                %dimension for all ray values
                %for gg=1:cone_num
                
                %preallocate empty vectors with limits for the  while loop
                max_loop=1+(SystemParam.maxbounce*SystemParam.scatnum);%max number of loops to run

                I_in=zeros(1,(max_loop));
                I_intrack=zeros(1,(max_loop));
                P_in=zeros((max_loop),2);
                Theta_in=zeros(1,(max_loop)); 
                Direction_in=zeros(1,(max_loop)); 
                %assign starting while loop values
                ray.Pow_enter(aa,xx,yy)=Ent_Pow(xx,yy,h);
                I_in(1,1)=Ent_Pow(xx,yy,h);
                Theta_in(1,1)=Theta(xx,yy,h);
                P_in(1,1:2)=[0,y0(xx,yy,h)];
                Direction_in(1,1)=1;
                %assign initial index values
                while_num=1;
                bounce_num=1;
                st_index=2;%starting index to allocate the first scattered ray we're tracking
                pow_check=0;
                MEAS0=sum(meas.inten);%initial measurement value
                while any(I_in) && while_num<=max_loop
            %%%%%%5%initial variable set up for while loop%%%%%
                    I_0=I_in(1,while_num);
                    I_intrack(1,while_num)=I_0;
                    pow_check=pow_check+I_0;
                    P_0=P_in(while_num,:);
                    Theta_0=Theta_in(1,while_num);
                    direction=Direction_in(1,while_num);
                    V_0=[direction*abs(cos(Theta_0)),sin(Theta_0)];
                    %initial measurement for use in later
                    %troubleshooting/difference checking
                    
                    meas_start0=sum(meas.inten);
                    b2h0=ray.b2hpow(aa,xx,yy);
                    cutoff0=ray.cutoffpow(aa,xx,yy);
                    trans0=ray.transmitted(aa,xx,yy);
                    abs0=ray.absorbed(aa,xx,yy);
                    back0=ray.backscat(aa,xx,yy);
                    smaabs0=ray.SMAabs(aa,xx,yy);
                    approx0=ray.approxpow(aa,xx,yy);
                    if I_0<SystemParam.I_min%check to make sure the while loop is worth it to keep running
                        ray.cutoffpow(aa,xx,yy)=ray.cutoffpow(aa,xx,yy)+I_0;%mark amount of light left after cutoff
                        I_in(1,while_num)=0;%redundant with the continue but just in case
                        while_num=while_num+1;%update the while number
                        continue%go to the next while loop
                    elseif length(I_0)>1 || length((P_0))>2 || length(Theta_0)>1 || length(V_0)>2
                        error('too big vector')
                    elseif isnan(I_0)
                        error('is NaN')                        
                    end

            %%%%%%traveling %%%%%%%%%%%%       
                    %traveling the length of the fiber
                    %start by recording the initial currently measured
                    %intensity
                    meas_start1=sum(meas.inten);
%                     %check transmission %should be in the travel function
                        [I_1,V_1,P_1,direction,meas,IT,Tally] = Traveling(I_0,V_0,P_0,direction,SystemParam,Bounds,Global_Index,Tally,meas);%temporary
                        %record various loss tallies
                        ray.SMAabs(aa,xx,yy)=IT.housi+ray.SMAabs(aa,xx,yy);
                        ray.absorbed(aa,xx,yy)=IT.absorbi+ray.absorbed(aa,xx,yy);
                        ray.backscat(aa,xx,yy)=IT.backi+ray.backscat(aa,xx,yy);
                        ray.cutoffpow(aa,xx,yy)=IT.cutoffi+ray.cutoffpow(aa,xx,yy);
                        ray.approxpow(aa,xx,yy)=IT.approxi+ray.approxpow(aa,xx,yy);
                        ray.b2hpow(aa,xx,yy)=IT.b2hi+ray.b2hpow(aa,xx,yy);
                        ray.transmitted(aa,xx,yy)=IT.transi+ray.transmitted(aa,xx,yy);
                     %after travel function
                    meas_dif=sum(meas.inten)-meas_start1;
                    Pow_out=[meas_dif,I_1,IT.backi,IT.absorbi,IT.housi,IT.cutoffi,IT.approxi];
                    [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_0,Pow_out,SystemParam,'Traveling',0,Global_Index,Tally);
                    ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
                    ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
            
                    %check if there's enough power to warrant going thru the end reflect
                    %function
                    if I_1<SystemParam.I_min
                        %store and record data. much more probably
                        ray.cutoffpow(aa,xx,yy)=ray.cutoffpow(aa,xx,yy)+I_1;%
                        %lots of stuff in here
                        %reset loop
                         I_in(1,while_num)=0;
                         while_num=while_num+1;
                         
                        continue%go to the next while loop
                    end
                    meas_dif_1=sum(meas.inten)-meas_start0;
                    b2h_1=ray.b2hpow(aa,xx,yy)-b2h0;
                    cutoff_1=ray.cutoffpow(aa,xx,yy)-cutoff0;
                    trans_1=ray.transmitted(aa,xx,yy)-trans0;
                    abs_1=ray.absorbed(aa,xx,yy)-abs0;
                    smaabs_1=ray.SMAabs(aa,xx,yy)-smaabs0;
                     approx_1=ray.approxpow(aa,xx,yy)-approx0;
                    back_1=ray.backscat(aa,xx,yy)-back0;
                    Pow_Out=[meas_dif_1,I_1,b2h_1,cutoff_1,trans_1,abs_1,smaabs_1,approx_1,back_1];
%                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
                    [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop step 0 to step 1',0,Global_Index,Tally);

                    
                    meas_start2=sum(meas.inten);%initial total measured before entering into the end reflect function
              %%%%%%%end reflect stuff%%%%%%%%%
                    if bounce_num<=SystemParam.maxbounce
                        %create the new rays for the loop, update the
                        %bounce number, take relevant measurements and
                        %tallies.
                        [I_scat,Theta_scat,I_trans,direction,b2h,bounce_num,Tally,meas]= EndReflect(I_1,V_1,P_1,SystemParam,direction,bounce_num,Global_Index,Tally,meas);  
                        %record and save
                        ray.transmitted(aa,xx,yy)=ray.transmitted(aa,xx,yy)+I_trans;
                        ray.b2hpow(aa,xx,yy)=b2h+ray.b2hpow(aa,xx,yy);%record and save
                        I_scat(I_scat<SystemParam.I_min)=0;%any of the scattered rays less than the minimum tracking power are set to 0
                        %create approprate indices for assigning the scattered values to the main while loop vectors
                        end_index=st_index+length(I_scat)-1;
                        indices_assign=st_index:1:end_index;
                        indices_assign=indices_assign(indices_assign<=max_loop);%make sure  the indices are under the max value
                        %sort the incoming vector to prioritize highest
                        %bvalues to be assigned to the vector
                        [I_sort_scat,sort_scat]=sort(I_scat,'descend');
                     
                        if ~isempty(indices_assign)%otherwise, nothing will need to be assigned to the loop vector
                            if end_index>max_loop%if there aren't enough spots to allocate all new I_sorted
                                highest_index=find(indices_assign==max_loop);%gives the total amount of spots available to allocate the sorted scI_Scat to
                                sort_scat=sort_scat(1:highest_index);
                            else%if theres plenty of space
                                highest_index=length(indices_assign);
                            end

                            %new loop assigned rays
                            I_in(1,indices_assign(1):indices_assign(highest_index))=I_scat(sort_scat);
                            Theta_in(1,indices_assign(1):indices_assign(highest_index))=Theta_scat(sort_scat);
                            P_in(indices_assign(1):indices_assign(highest_index),1)=P_1(1);P_in(indices_assign(1):indices_assign(end),2)=P_1(2);%all will have same initial point
                            Direction_in(1,indices_assign(1):indices_assign(highest_index))=direction;
                            sum_forwardscat=sum(I_scat(sort_scat));
                        else
                            sum_forwardscat=0;

                        end
                        %update the starting index
                        st_index=end_index+1;
                    else%if we've already had the max number of bounces, just record the transmitted and/or lost light
                        [I_scat,~,I_trans,~,b2h,~,Tally,meas]= EndReflect(I_1,V_1,P_1,SystemParam,direction,bounce_num,Global_Index,Tally,meas);
                        ray.transmitted(aa,xx,yy)=ray.transmitted(aa,xx,yy)+I_trans;
                        ray.b2hpow(aa,xx,yy)=b2h+ray.b2hpow(aa,xx,yy);%record and save
                        %record what isn't being tracked
                        ray.cutoffpow(aa,xx,yy)=sum(I_scat)+ray.cutoffpow(aa,xx,yy);
                        sum_forwardscat=0;
                    end
                    %check if there's any major differences in the total
                    meas_dif=sum(meas.inten)-meas_start2;
                    [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_1,[meas_dif,sum(I_scat),I_trans,b2h],SystemParam,'End Reflect',0,Global_Index,Tally);
                    ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
                    ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
                    
                    %set the just used I_in to 0
                    I_in(1,while_num)=0;

                    %update the while counter
                    while_num=while_num+1;
                    if while_num<(max_loop)%more than one while loop remaining
                        %sort all of the remaining in descending order so the
                        %highest intensity ones are prioritized
                        [I_in(1,while_num:end),sort_Index]=sort(I_in(while_num:end),'descend');
                        sort_Index=sort_Index+(while_num-1);
                        Theta_in(1,while_num:end)=Theta_in(1,sort_Index);
                        P_in(while_num:end,1)=P_in(sort_Index,1);P_in(while_num:end,2)=P_in(sort_Index,2);
                        Direction_in(1,while_num:end)=Direction_in(1,sort_Index);
                    end
                    
                    %if the loop cuts off with anything remaining in it
                    if while_num==(2+(SystemParam.maxbounce*SystemParam.scatnum)) && any(I_in)%these conditions shouldn't occur BUT
                        ray.cutoffpow(aa,xx,yy)=sum(I_in)+ray.cutoffpow(aa,xx,yy);
                    end

                    meas_dif_final=sum(meas.inten)-meas_start0;
                    b2h_final=ray.b2hpow(aa,xx,yy)-b2h0;
                    cutoff_final=ray.cutoffpow(aa,xx,yy)-cutoff0;
                    trans_final=ray.transmitted(aa,xx,yy)-trans0;
                    abs_final=ray.absorbed(aa,xx,yy)-abs0;
                    smaabs_final=ray.SMAabs(aa,xx,yy)-smaabs0;
                    approx_final=ray.approxpow(aa,xx,yy)-approx0;
                    back_final=ray.backscat(aa,xx,yy)-back0;
%                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
                    Pow_out=[sum_forwardscat,meas_dif_final,b2h_final,cutoff_final,trans_final,abs_final,smaabs_final,approx_final,back_final];
                    [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop',0,Global_Index,Tally);

                end
                %summing all of the measured side emitted power
                
                ray.pow_side(aa,xx,yy)=sum(meas.inten)-MEAS0;
                total_pow_use=[ray.pow_side(aa,xx,yy),ray.transmitted(aa,xx,yy),ray.SMAabs(aa,xx,yy),ray.b2hpow(aa,xx,yy),ray.backscat(aa,xx,yy),ray.absorbed(aa,xx,yy),ray.cutoffpow(aa,xx,yy)];
                
                %tracking differences
                [Tally,Differenceamount,Diffamountpos] =  DifTrack(ray.Pow_enter(aa,xx,yy),total_pow_use,SystemParam,'total while loop',0,Global_Index,Tally);
                ray.remaininglosses(aa,xx,yy)=Differenceamount;
                %store measured data
                %only include values that exist
                indexes=find(meas.inten);
                 ray.meas_inten(aa,xx,yy)={meas.inten(indexes)};
                 ray.meas_points(aa,xx,yy)={meas.points(indexes,:)};
                %place to assign summed values to (xx*yy) sized column
                %vector
                prev_xxnum=xx-1;
                xxyy_ind=yy+(prev_xxnum*b);
                aa_Rand(1).meas_inten(xxyy_ind,aa)={meas.inten(indexes)};%assign measured intensity to a row vector
                aa_Rand(1).meas_points(xxyy_ind,aa)=ray.meas_points(aa,xx,yy);%assign measured points to a row vector

 
                end
                
            end
            
            
            
            meas_point_aa=cell2mat(aa_Rand(1).meas_points(:,aa));
            meas_inten_aa=cell2mat(aa_Rand(1).meas_inten(:,aa));
                if any(meas_inten_aa)
                    disp('any')
                       [XVEC,Y, lengthPlot,pow_side_total] = Bins(meas_point_aa, meas_inten_aa, SystemParam.division, xlen, r_fiber,SystemParam);
                    use_index=find(meas.points(:,1)>=(2.5*10^4));%indexes of all of the measurement points that correspond to measurable light after the sma
                    water_st_index=find(meas_point_aa(:,1)>=SystemParam.waterstart);
                    SMA_index=find(meas_point_aa(:,1)<=SystemParam.SMA_totallength);
                    pow_side_use=sum(meas_inten_aa(use_index));
                    SMA_pow_side=sum(meas_inten_aa(SMA_index));
                    pow_side_waterst=sum(meas_inten_aa(water_st_index));
                    %performance metric parameters
                    aa_Rand(1).Y(aa,1)=Y;
                    Ylocal=cell2mat(Y);%local version of the I(x) vector
                    aa_Rand(1).UC(1,aa)=Ylocal(i6)/Ylocal(i1);
                else
                pow_side_total=0;
                pow_side_use=0;
                SMA_pow_side=0;
                pow_side_waterst=0;
                %performance metric parameters
                Ylocal=zeros(1,(floor((SystemParam.xlen/10e3)/(SystemParam.division/10e3))+1));%local version of the I(x) vector
                aa_Rand(1).UC(1,aa)=0;
                aa_Rand(1).Y(aa,1)={Ylocal};
                end
        %update aa_lim storage
        %summing storage
        
        aa_Rand(1).Pow_enter(1,aa)=sum(ray.Pow_enter(aa,:,:),'all');
        aa_Rand(1).transmitted(1,aa)=sum(ray.transmitted(aa,:,:),'all');
        aa_Rand(1).RatioIsIt(aa,xx,yy)=aa_Rand(1).pow_side(1,aa)/aa_Rand(1).transmitted(1,aa);
        aa_Rand(1).pow_side(1,aa)=pow_side_total;
        aa_Rand(1).pow_side_use(1,aa)=pow_side_use;
        aa_Rand(1).SMA_pow_side(1,aa)=SMA_pow_side;
        aa_Rand(1).pow_side_waterst(1,aa)=pow_side_waterst(1,aa);
        aa_Rand(1).absorbed(1,aa)=sum(ray.absorbed(aa,:,:),'all');
        aa_Rand(1).backscat(1,aa)=sum(ray.backscat(aa,:,:),'all');
        aa_Rand(1).SMAabs(1,aa)=sum(ray.SMAabs(aa,:,:),'all');
        aa_Rand(1).approxpow(1,aa)=sum(ray.approxpow(aa,:,:),'all');
        aa_Rand(1).approxpow_dif(1,aa)=sum(ray.approxpow_dif(aa,:,:),'all');
        aa_Rand(1).approxpow_pos(1,aa)=sum(ray.approxpow_pos(aa,:,:),'all');
        aa_Rand(1).b2hpow(1,aa)=sum(ray.b2hpow(aa,:,:),'all');
        aa_Rand(1).cutoffpow(1,aa)=sum(ray.cutoffpow(aa,:,:),'all');
        %average metrics
        aa_Rand(1).RatioIsIt(1,aa)= aa_Rand(1).pow_side_use(1,aa)./aa_Rand(1).transmitted(1,aa);

        
        %aa_Rand(2).Y(aa,1)={std(cell2mat(ray.Y(aa,:,:)))};
%             
        %actually have to calculate this one
        aa_Rand(1).remaininglosses=aa_Rand(1).Pow_enter(1,aa)-(aa_Rand(1).pow_side(1,aa)+aa_Rand(1).transmitted(1,aa)+aa_Rand(1).SMAabs(1,aa)+aa_Rand(1).b2hpow(1,aa)+aa_Rand(1).backscat+ aa_Rand(1).absorbed(1,aa)+aa_Rand(1).cutoffpow(1,aa));
% 
%             %include a place to plot instant aa_lim
            
        end
        %storing average data 
        FibIt(1).Pow_enter(iteration,h)=mean(aa_Rand(1).Pow_enter(1,:));
        FibIt(1).transmitted(iteration,h)=mean(aa_Rand(1).transmitted(1,:));
        FibIt(1).pow_side(1,aa)=mean(aa_Rand(1).pow_side(1,:));
        FibIt(1).pow_side_use(1,aa)=mean(aa_Rand(1).pow_side_use(1,:));
        FibIt(1).SMA_pow_side(1,aa)=mean(aa_Rand(1).SMA_pow_side(1,:));
        FibIt(1).pow_side_waterst(1,aa)=mean(aa_Rand(1).pow_side_waterst(1,:));
        FibIt(1).absorbed(1,aa)=mean(aa_Rand(1).absorbed(1,:));
        FibIt(1).backscat(1,aa)=mean(aa_Rand(1).backscat(1,:));
        FibIt(1).SMAabs(1,aa)=mean(aa_Rand(1).SMAabs(1,:));
        FibIt(1).approxpow(1,aa)=mean(aa_Rand(1).approxpow(1,:));
        FibIt(1).approxpow_dif(1,aa)=mean(aa_Rand(1).approxpow_dif(1,:));
        FibIt(1).approxpow_pos(1,aa)=mean(aa_Rand(1).approxpow_pos(1,:));
        FibIt(1).b2hpow(1,aa)=mean(aa_Rand(1).b2hpow(1,:));
        FibIt(1).cutoffpow(1,aa)=mean(aa_Rand(1).cutoffpow(1,:));
        FibIt(1).remaininglosses=mean(aa_Rand(1).remaininglosses);

        %std deviations
        FibIt(2).Pow_enter(iteration,h)=std(aa_Rand(1).Pow_enter(1,:));
        FibIt(2).transmitted(iteration,h)=std(aa_Rand(1).transmitted(1,:));
        FibIt(2).pow_side(iteration,h)=std(aa_Rand(1).pow_side(1,:));
        FibIt(2).pow_side_use(iteration,h)=std(aa_Rand(1).pow_side_use(1,:));
        FibIt(2).SMA_pow_side(iteration,h)=std(aa_Rand(1).SMA_pow_side(1,:));
        FibIt(2).pow_side_waterst(iteration,h)=std(aa_Rand(1).pow_side_waterst(1,:));
        FibIt(2).absorbed(iteration,h)=std(aa_Rand(1).absorbed(1,:));
        FibIt(2).backscat(iteration,h)=std(aa_Rand(1).backscat(1,:));
        FibIt(2).SMAabs(iteration,h)=std(aa_Rand(1).SMAabs(1,:));
        FibIt(2).approxpow(iteration,h)=std(aa_Rand(1).approxpow(1,:));
        FibIt(2).approxpow_dif(iteration,h)=std(aa_Rand(1).approxpow_dif(1,:));
        FibIt(2).approxpow_pos(iteration,h)=std(aa_Rand(1).approxpow_pos(1,:));
        FibIt(2).b2hpow(iteration,h)=std(aa_Rand(1).b2hpow(1,:));
        FibIt(2).cutoffpow(iteration,h)=std(aa_Rand(1).cutoffpow(1,:));
        FibIt(2).remaininglosses=std(aa_Rand(1).remaininglosses);

        %average metrics
        if aa_lim==1
             FibIt(1).Y(iteration,h)=aa_Rand(1).Y;
        else
            FibIt(1).Y(iteration,h)={mean(cell2mat(aa_Rand(1).Y),2)};
            FibIt(1).Y(iteration,h)={mean(cell2mat(aa_Rand(1).Y),2)};
        end


        FibIt(1).UC(iteration,h)=mean(aa_Rand(1).UC(1,:));
        FibIt(1).RatioIsIt(iteration,h)=mean(aa_Rand(1).RatioIsIt(1,:));
        FibIt(2).UC(iteration,h)=std(aa_Rand(1).UC(1,:));
        FibIt(2).RatioIsIt(iteration,h)=std(aa_Rand(1).RatioIsIt(1,:));
                

         %storing the larger data for maybe checking the trends out with
         SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
         SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
         SumVal.transmitted(iteration,h)={ray.transmitted};
         SumVal.pow_side(iteration,h)={ray.pow_side};
         SumVal.pow_side_use(iteration,h)={ray.pow_side_use};
         SumVal.SMA_pow_side(iteration,h)={ray.SMA_pow_side};
         SumVal.pow_side_waterst(iteration,h)={ray.pow_side_waterst};
         SumVal.absorbed(iteration,h)={ray.absorbed};
         SumVal.backscat(iteration,h)={ray.backscat};
         SumVal.SMAabs(iteration,h)={ray.SMAabs};
         SumVal.approxpow_dif(iteration,h)={ray.approxpow_dif};%sum of power differences due to difference in approximations
         SumVal.approxpow_pos(iteration,h)={ray.approxpow_pos}; %absolute value of all of the power difference approximationsSumVal.b2hpow(iteration,h)={ray.;
         SumVal.cutoffpow(iteration,h)={ray.cutoffpow};
         SumVal.remaininglosses(iteration,h)={aa_Rand(1).remaininglosses};
         SumVal.UC(iteration,h)={ray.UC};
         SumVal.RatioIsIt(iteration,h)={ray.RatioIsIt};
         SumVal.Y(iteration,h)={ray.Y};
         
         %plotting the pie chart tracking stuff
         %percent of each thing for power balance
figure(1)
Power_Track=[aa_Rand(1).pow_side(1,aa),aa_Rand(1).transmitted(1,aa),aa_Rand(1).SMAabs(1,aa),aa_Rand(1).b2hpow(1,aa),aa_Rand(1).backscat(1,aa),aa_Rand(1).absorbed(1,aa),aa_Rand(1).cutoffpow(1,1),aa_Rand(1).remaininglosses(1,1)];
Power_Track_lab=["pow_side" "transmitted" "SMAabs" "b2hpow" "backscat" "absorbed" "cutoffpow" "remaininglosses"];
pie(Power_Track,Power_Track_lab)

figure(2)
Power_Track_comp=[aa_Rand(1).pow_side_use(1,aa),aa_Rand(1).SMA_pow_side(1,aa),aa_Rand(1).transmitted(1,1),aa_Rand(1).SMAabs(1,1),aa_Rand(1).b2hpow(1,1),aa_Rand(1).backscat(1,1),aa_Rand(1).absorbed(1,1),aa_Rand(1).cutoffpow(1,1),aa_Rand(1).remaininglosses(1,1)];
Power_Track_comp_lab=["pow_side_use" "SMA_pow_side" "transmitted" "SMAabs" "b2hpow" "backscat" "absorbed" "cutoffpow" "remaininglosses"];
pie(Power_Track_comp,Power_Track_comp_lab)
% 
%plotting the I(x) distribution data for each iteration
figure(3)
%title and legend defined earlier lines 95 is through 140ish
Y_current=cell2mat(FibIt(1).Y(iteration,h));%current I(x) within this iteration and fiber
plot(XVEC,Y_current)
title(Title_Main) %defined earlier 
hold on
if yes_itname==1%if we need a legend bc there's multiple iterations
    legend(legend_Main)
end



%%%%% this is all just extracting the data from the Tally struct to find where the biggest differences in coming from for trouble shooting
%%%% i didn't comment this well nor is this probably the best method to
%%%% extract the data. but it works for now
%figuring out where the gaps are, sort the strings of the difference values
[m1,m2,m3,m4]=size((Tally.big_dif_fun_main));%[iteration, fiber number, aa limit, number of times thing got assigned to the tally]
[f1,f2,f3,f4,f5,f6]=size((Tally.big_dif_fun_func));%[iteration, fiber number, aa limit,xx,yy, number of times thing got assigned to the tally]

maincount=numel(Tally.big_dif_count_meas_main);%total number of times we count a "big difference" in the main code
funcount=numel(Tally.big_dif_count_meas_func);%total number of times we count a "big difference" within each function
IO_maincount=numel(Tally.IO_count_meas_main);
IO_funcount=numel(Tally.IO_count_meas_func);

bigmain=strings(1,maincount);%empty vector for identifying the part of the code that has the big difference
dif_main=zeros(1,maincount);%empty vector for recording the big difference within the main
dif_fun=zeros(1,funcount);%empty vector for identifying the function that has the big difference
bigfun=strings(1,funcount);%empty vector for recording the big difference within the function
cmain=1;%index
cfun=1;%index for the main code
%for the in and out check
IO_bigmain=strings(1,IO_maincount);%empty vector for identifying the part of the code that has the big difference
IO_dif_main=zeros(1,IO_maincount);%empty vector for recording the big difference within the main
IO_dif_fun=zeros(1,IO_funcount);%empty vector for identifying the function that has the big difference
IO_bigfun=strings(1,IO_funcount);%empty vector for recording the big difference within the function
IO_cmain=1;%index
IO_cfun=1;%index for the main code

for i=1:m1 %iteration number
    for j=1:m2 %fiber number
        for k=1:m3 %aa number
            %set up emtry string vector
           count_inst=Tally.big_dif_count_main(i,j,k);%find maximum index where there is a location in main that had
           IO_count_inst=Tally.IO_count_main(i,j,k);%find maximum index where there is a location in main that had

           %a big difference
            for l=1:(count_inst)
                if ischar(Tally.big_dif_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
                   
                    bigmain(1,cmain)=Tally.big_dif_fun_main{i,j,k,l}{1,1};%assign the string value
                    dif_main(1,cmain)=cell2mat(Tally.big_dif_count_meas_main{i,j,k,l}); %assign the difference value
                    cmain=cmain+1; %update the index
                end
            end
            for l=1:(IO_count_inst)
                if ischar(Tally.IO_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
                   
                    IO_bigmain(1,IO_cmain)=Tally.IO_fun_main{i,j,k,l}{1,1};%assign the string value
                    IO_dif_main(1,IO_cmain)=cell2mat(Tally.IO_count_meas_main{i,j,k,l}); %assign the difference value
                    IO_cmain=IO_cmain+1; %update the index
                end
            end
            %assign
            for l=1:f4
                for m=1:f5
                    %find each location where there is a function that had
                    %a big difference
                    count_inst=find(Tally.big_dif_fun_func{i,j,k,l,m});%
                    IO_count_inst=Tally.IO_count_func(i,j,k,l,m);

                    for o=1:length(count_inst)
                        n=count_inst(o);%the index
                            %double check the difference function isn't
                            %empty
                            if isempty(Tally.big_dif_fun_func{i,j,k,l,m,n})
                                fun_str_inst="";
                            else
                                fun_str_inst=Tally.big_dif_fun_func{i,j,k,l,m,n};
                                 dif_fun(1,cfun)=Tally.big_dif_count_meas_func{i,j,k,l,m,n};
                                 bigfun(1,cfun)=fun_str_inst;
                                 cfun=cfun+1;

                            end
                    end
                    if IO_count_inst~=0
                        for o=1:length(IO_count_inst)
                            n=IO_count_inst(o);%the index
                                %double check the difference function isn't
                                %empty
                                if isempty(Tally.IO_fun_func{i,j,k,l,m,n})
                                    fun_str_inst="";
                                else
                                    fun_str_inst=Tally.IO_fun_func{i,j,k,l,m,n};
                                    IO_dif_fun(1,IO_cfun)=Tally.IO_count_meas_func{i,j,k,l,m,n};
                                    IO_bigfun(1,IO_cfun)=fun_str_inst;
                                    IO_cfun=IO_cfun+1;

                                end
                        end
                    end
                end
            end
        end
    end
end
idfun=find(bigfun~="");
    if sum(dif_fun,'all')==sum(dif_fun(idfun),'all')
        bigfun=bigfun(idfun);
        fun_cats=unique(bigfun);%create a category of each unique function reported with a big dif
        trouble_functions=categorical(bigfun,fun_cats);%,fun_categorical);
    else
        error('figure out the indexing issues')
    end
main_cats=unique(bigmain);%create a category of each unique function reported with a big dif
summary(trouble_functions)%get a summary of all of the trouble functions
trouble_main=categorical(bigmain,main_cats);%main_categorical);
summary(trouble_main) %get a summary of all of the trouble areas in the main code

IO_idfun=find(IO_bigfun~="");
    if sum(IO_dif_fun,'all')==sum(IO_dif_fun(IO_idfun),'all')
        IO_bigfun=IO_bigfun(IO_idfun);
        IO_fun_cats=unique(IO_bigfun);%create a category of each unique function reported with a big dif
        IO_trouble_functions=categorical(IO_bigfun,IO_fun_cats);%,fun_categorical);
    else
        error('figure out the indexing issues')
    end
IO_idmain=find(IO_bigmain~="");    
IO_main_cats=unique(IO_bigmain(IO_idmain));%create a category of each unique function reported with a big dif
summary(IO_trouble_functions)%get a summary of all of the trouble functions
IO_trouble_main=categorical(IO_bigmain,IO_main_cats);%main_categorical);
summary(IO_trouble_main) %get a summary of all of the trouble areas in the main code

%create a table with summary values of  issues caused by each
%function
%empty cell for all of the categories of each big_dif item
dif_f=cell(length(fun_cats),1);
dif_m=cell(length(main_cats),1);
IO_dif_f=cell(length(fun_cats),1);
IO_dif_m=cell(length(main_cats),1);
%empty vectory for all of the categories of each big_dif item
dif_contrib_m=zeros(length(main_cats),2);
dif_contrib_f=zeros(length(fun_cats),2);

IO_dif_contrib_m=zeros(length(IO_main_cats),2);
IO_dif_contrib_f=zeros(length(IO_fun_cats),2);
for i=1:length(main_cats)%for each categorical in the main function
    dif_m(1,i)={dif_main(1,(trouble_main==main_cats(i)))};% max a cell 
    dif_contrib_m(i,1)=sum(cell2mat(dif_m(1,i)),'all');
    dif_contrib_m(i,2)=sum(abs(cell2mat(dif_m(1,i))),'all');
end
for i=1:length(fun_cats)
    dif_f(1,i)={dif_fun(1,(trouble_functions==fun_cats(i)))};
    dif_contrib_f(i,1)=sum(cell2mat(dif_f(1,i)),'all');
    dif_contrib_f(i,2)=sum(abs(cell2mat(dif_f(1,i))),'all');
end

for i=1:length(IO_main_cats)%for each categorical in the main function
    IO_dif_m(1,i)={IO_dif_main(1,(IO_trouble_main==IO_main_cats(i)))};% max a cell 
    IO_dif_contrib_m(i,1)=sum(cell2mat(IO_dif_m(1,i)),'all');
    IO_dif_contrib_m(i,2)=sum(abs(cell2mat(IO_dif_m(1,i))),'all');
end
for i=1:length(IO_fun_cats)
    IO_dif_f(1,i)={IO_dif_fun(1,(IO_trouble_functions==IO_fun_cats(i)))};
    IO_dif_contrib_f(i,1)=sum(cell2mat(IO_dif_f(1,i)),'all');
    IO_dif_contrib_f(i,2)=sum(abs(cell2mat(IO_dif_f(1,i))),'all');
end


Changed=dif_contrib_m(:,1);%difference value vector for the table
Total=dif_contrib_m(:,2);%absolute value vector for the table
Function=transpose(main_cats); %function
TabM=table(Changed,Total,'RowNames',Function);
Changed=dif_contrib_f(:,1);%difference value vector for the table
Total=dif_contrib_f(:,2);%absolute value vector for the table
Function=transpose(fun_cats);
TabF=table(Changed,Total,'RowNames',Function);


IO_Changed=IO_dif_contrib_m(:,1);%difference value vector for the table
IO_Total=IO_dif_contrib_m(:,2);%absolute value vector for the table
IO_Function=transpose(IO_main_cats); %function
IO_TabM=table(IO_Changed,IO_Total,'RowNames',IO_Function);
IO_Changed=IO_dif_contrib_f(:,1);%difference value vector for the table
IO_Total=IO_dif_contrib_f(:,2);%absolute value vector for the table
IO_Function=transpose(IO_fun_cats);
IO_TabF=table(IO_Changed,IO_Total,'RowNames',IO_Function);


%%%%%%writing the data to an excel file
    %%%%SAVING EACH ITERATION, each fiber DATA
    %saving
        if (writeToFile == 1)
            %documentation of variables used in the 1st row
            if c>1 %if multiple fibers
                description=description + " fiber no. " +h;
            end
            writematrix(description, filename,'Sheet', sheet_num + (h-1),'Range', "A1");
            writematrix("Separation Distance: " + SystemParam.led_distance + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "B1");
            writematrix("Wavelength: " + SystemParam.uv_wavelength + "(nm)", filename,'Sheet', sheet_num + (h-1),'Range', "C1");
            writematrix("Angle resolution: " + SystemParam.angnum, filename,'Sheet', sheet_num + (h-1),'Range', "E1");
            writematrix("initial intensity Entering the Fiber: " + SystemParam.I_init + " (mW)", filename,'Sheet', sheet_num + (h-1),'Range', "F1");
    %         writematrix("Polymer coat: " + SystemParam.poly_coat, filename,'Sheet', sheet_num + (h-1),'Range', "G1");
    %         writematrix("Polymer thickness (if true): " + polymer_thickness + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "H"+ num2str(iteration));
    %         writematrix("Nanoparticle Density: " + np_density, filename,'Sheet', sheet_num + (h-1),'Range', "I1");
    %         writematrix("Nanoparticle Layer Thickness: " + np_thickness + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "J1");
    %         writematrix("Nanoparticle Diameter: " + SystemParam.particle_dmtr + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "K1");
            %writematrix("# of Mie vectors considered: " + num_vectors, filename,'Sheet', sheet_num + (h-1),'Range', "L1");   
            %writematrix("Measurement distance: " + SystemParam.meas_distance + ( "um"), filename,'Sheet', sheet_num + (h-1),'Range', "M1");
            writematrix(ext_media, filename,'Sheet', sheet_num + (h-1),'Range', "N1");
            writematrix("# of for loops used to smooth random noise: " + aa_lim, filename,'Sheet', sheet_num + (h-1),'Range', "O1");
            writematrix("Length of fiber division in plot: " + (SystemParam.division/10e3) + " (cm)", filename,'Sheet', sheet_num + (h-1),'Range', "P1");
            writematrix("Photons Considered: " +SystemParam.raynum + ", Angles Considered: "+  SystemParam.angnum, filename,'Sheet', sheet_num + (h-1),'Range', "Q1");
            writematrix("Number Fibers: " +SystemParam.num_fib, filename,'Sheet', sheet_num + (h-1),'Range', "R1");
            housepow=1;
            %iteration values
            writematrix("Fiber length: " + (SystemParam.xlen/10000) + " (cm)", filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)-3));
            writematrix("Radius of fiber: " + SystemParam.r_fiber + " (um)", filename,'Sheet', sheet_num + (h-1),'Range', "B" + num2str((6*iteration)-3));
            writematrix("Distance of LED from Fiber: " +LED_d + " (mm)", filename,'Sheet', sheet_num + (h-1),'Range', "C"+ num2str((6*iteration)-3));
            %documentation of the data generated in the iteration
            writematrix("Total light entering the fiber: " + FibIt(1).Pow_enter(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "E" + num2str((6*iteration)-3));
            writematrix("Total transmitted light: " + FibIt(1).transmitted(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "F" + num2str((6*iteration)-3));
            writematrix("Total side emitted light: " + FibIt(1).pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "G" + num2str((6*iteration)-3));
            writematrix("Total useable side emitted light: " + FibIt(1).pow_side_use(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "H" + num2str((6*iteration)-3));
            writematrix("Total side emitted light w/in the SMA: " + FibIt(1).SMA_pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "I" + num2str((6*iteration)-3));
            writematrix("Total adsorbed light: " + FibIt(1).absorbed(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "J" + num2str((6*iteration)-3));
            writematrix("Total back scattered light: " + FibIt(1).backscat(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "K" + num2str((6*iteration)-3));
            writematrix("Total SMA adsorbed light: " + FibIt(1).SMAabs(iteration,h)+ " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "L" + num2str((6*iteration)-3));
            writematrix("Total approximation loss light: " + FibIt(1).approxpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "M" + num2str((6*iteration)-3));
            writematrix("Total back 2housing loss light: " + FibIt(1).b2hpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "N" + num2str((6*iteration)-3));
            writematrix("Total cutoff loss light: " + FibIt(1).cutoffpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "O" + num2str((6*iteration)-3));
            writematrix("Total remaining light loss: " + FibIt(1).remaininglosses(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "P" + num2str((6*iteration)-3));
            writematrix("Avg Uniformity Cofficient: " + FibIt(1).UC(iteration,h), filename,'Sheet', sheet_num + (h-1),'Range', "Q" + num2str((6*iteration)-3));
            writematrix("Avg Is/It Ratio: " + FibIt(1).RatioIsIt(iteration,h) , filename,'Sheet', sheet_num + (h-1),'Range', "R" + num2str((6*iteration)-3));
            %documentation of the std deviation of the generated data
            writematrix(FibIt(2).Pow_enter(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "E" + num2str((6*iteration)-2));
            writematrix(FibIt(2).transmitted(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "F" + num2str((6*iteration)-2));
            writematrix(FibIt(2).pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "G" + num2str((6*iteration)-2));
            writematrix(FibIt(2).pow_side_use(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "H" + num2str((6*iteration)-2));
            writematrix(FibIt(2).SMA_pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "I" + num2str((6*iteration)-2));
            writematrix(FibIt(2).absorbed(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "J" + num2str((6*iteration)-2));
            writematrix(FibIt(2).backscat(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "K" + num2str((6*iteration)-2));
            writematrix(FibIt(2).SMAabs(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "L" + num2str((6*iteration)-2));
            writematrix(FibIt(2).approxpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "M" + num2str((6*iteration)-2));
            writematrix(FibIt(2).b2hpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "N" + num2str((6*iteration)-2));
            writematrix(FibIt(2).cutoffpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "O" + num2str((6*iteration)-2));
            writematrix(FibIt(2).remaininglosses(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "P" + num2str((6*iteration)-2));
            writematrix(FibIt(2).UC(iteration,h), filename,'Sheet', sheet_num + (h-1),'Range', "Q" + num2str((6*iteration)-2));
            writematrix(FibIt(2).RatioIsIt(iteration,h) , filename,'Sheet', sheet_num + (h-1),'Range', "R" + num2str((6*iteration)-2));
            %side emission over a the fiber distance and iteration
            writematrix(XVEC,filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)-1));
            writematrix(cell2mat(FibIt(1).Y(iteration,h)),filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)));
            writematrix(cell2mat(FibIt(2).Y(iteration,h)),filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)+1));

        end
    end%%%%%%%%%%

end
hold off


