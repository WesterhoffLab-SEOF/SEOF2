%%simplest possible fiber model%%
%changes frp, prev are just an update to the binning model to be more
%representative of the x data 

close all; clear all;%fresh slate
addpath('Utilities');%functions go into utilities folder
addpath('Data Structs');%structs go into this folder

%%%%%%%%%%%%SETTING UP THE SIMULATION%%%%%%%%%%%%%%
%systemParam, direct input, can change these
SystemParam = setupSystemParams();
%data display variables... what do you want to see generated?
writeToFile =1;   % set to 1 to save generated results to excel file, set to 0 otherwise
%filename = 'C:\Users\ndshapi1.ASURITE\OneDrive - Arizona State University\SEOF_Raytracing 082823\data\medswings50cm041525_air.xlsx';%%C:\Users\Nora Shapiro\OneDrive - Arizona State University\SEOF_Raytracing 082823\data\model_50cm_111023.xlsx'; %name of file %MUST INCLUDE ADDRESS OF DESTINATION FOLDER    
%Make sure the filename is in a writable location and does not require
%admin access
filename='C:\Users\etwes\OneDrive\Documents\test.xlsx';%
sheet_num=1;%Excel Sheet number t
description="21x21 rays, 100mW 170deg, 70 deg half angle, led emission, 1mm sep dist. 50.5* water 4.5cm with unfilled sma connector*. nfibi [2.5,3.75,5]iE-07,rcoef[0.7,0.725,.75,0.775,0.8],ang=pi/[15,18],nmetal[1.3,1.325,1.35],scat=[11,15]";%description of the changes made in this file



%%%SETTING UP PARAMETERS NEEDED FOR ITERATION AND DATA DISPLAY%%%%
%all of these variables can be turned into a vector of different values so
%it will automatically re-run the model with a new variable.
%CAN ONLY ITERATE ON ONE AT A TIME FOR NOW
rfiber=SystemParam.r_fiber;%[125,200,250,500];%;%SystemParam.r_fiber;%[1500/2,1500];%[125,200,250,500];%[125,200,250,500];%[125,250,500]%um alternate
x_len=SystemParam.xlen;%[12.5,50.5].*10^4;%([10,15,20,25,30,35,40,45,50]+2)*10^4;%cm->[um] alternate
xlen=SystemParam.xlen;
led_dist=SystemParam.led_distance;%[0,1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000];% mm alternate[0,1,2,3,4,5];% mm alternate
nfib=1.5+1i;%SystemParam.n1;%1.5+([1E-07,1.2E-07,1.5E-07,2.5E-07,5E-07].*1i);%description of the changes made in this file
%1.5+([0.5E-07,1E-07,1.2E-07,1.3E-07,1.4E-07,1.5E-07,1.6E-07,1.7E-07,1.8E-07,2E-07,2.5E-07,5E-07].*1i);%SystemParam.n1;%1.5+(1i.*logspace(-7.301,-4.6021,15));%1.5+(1.5E-07*1i);%%1.5+[1.5E-07,2.5E-07,5E-07].*1i;%1.5+[1.5E-07].*1i;%(linspace(1E-08,1E-07,4).*1i);%1.5+(1i.*1.06*10^-7);%linspace(6*10^-8,3*10^-7,10);%SystemParam.n1;
% nmetal_real=1.3+linspace(-0.3,0.1,5);
% nmetal_comp=(2.3+linspace(-0.2,0.2,5))*1i;
% [nmetal_r,nmetal_i]=meshgrid(nmetal_real,nmetal_comp);
n_metal=1.325+2.3i;%SystemParam.n_metal;%[1.3,1.35,1.37,1.41,1.42,1.45,1.5]+2.3i;%;[1.38,1.39,1.4]+2.3i;%[(1.39+2.3i),(1.365+0.01i),(0.22+3.2i)];%%[1.3,1.45,1.475,1.5]+2.3i;%%(1.4,1.45,1.46,1.47):0.02:1.48)+(1i*2.3);%;%[(1.46+2.3i),(1.365+0.01i),0.22+3.2i]%(1.445:0.005:1.465)+2.3i;%%[1.3,1.45,1.475,1.5]+2.3i;%%(1.32:0.02:1.48)+(1i*2.3);%nmetal_r+nmetal_i;
ray_sqrt=SystemParam.angnum;%[21,23,25,27,29,31,41,51,101];%
scatter_num=7;%SystemParam.scatnum;%[11,15,19];%SystemParam.scatnum;%[3,7,9,11,13,21,27];%SystemParam.scatnum;%[3,7,9,11,13,21,27];%%13:2:27;%
dx_num=SystemParam.cont_dx;%[100000,10000,1000,100,10,1];%
I_min_var=SystemParam.I_min;%SystemParam.photon_min*10.^([15,12,9,6,3,1]);%
sma_num =SystemParam.housing_bounce;%[1,2,3,4,5];%
%ang_vec=linspace(pi/7.2,pi/6.8,9);
max_scatter_ang=pi/15;%SystemParam.scat_ang_max;%[pi/3,pi/7,pi/18];%SystemParam.scat_ang_max;%[pi/2.1,pi/6.75,pi/5,pi/8,pi/18];%[pi/6.75,pi/6.9,pi/7,pi/7.1,pi/7.25];%%linspace(pi/7.2,pi/5.8,7);%[pi/7.2,pi/6.8,pi/6.4];%linspace(pi/6,pi/2,5);%[pi/18,pi/15,pi/12,pi/9,pi/6,pi/4,pi/3,pi*2/5,pi*3/7,pi/2];% [pi/18,pi/6,pi/4,pi/3,pi/2];%
bounce_var=SystemParam.maxbounce;%[1,2,3,4,5];%
raysc_coeff=0.7;%SystemParam.r_coeff;%[0.9,0.6,0.3];%[0.01,0.5,0.9];%[0.1:0.025:0.2,0.3];%SystemParam.r_coeff;%[0.1,0.125,0.15,0.175,0.2];%[0.1,0.25,0.5];%%coefficient for what percentage of the loss is attributed to scattering vs absorption
sma_fill_length=SystemParam.SMA_fill;%(1.25:0.25:2.25).*10^4;%SystemParam.SMA_fill;%(1.25:.25:2.5).*10^4;%SystemParam.SMA_fill;
water_status=SystemParam.waterInterface;%[1,0];%[1,0];%
%include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vector of the size of each possible iterateable vector
%lengths of each iterable variable
lengths_it=[length(rfiber),length(x_len),length(led_dist),length(nfib),length(n_metal),length(ray_sqrt),length(scatter_num),...
    length(dx_num),length(I_min_var),length(sma_num),length(max_scatter_ang),length(bounce_var),length(raysc_coeff),length(sma_fill_length),length(water_status)];%include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names_it=['rfiber','xlen','led_dist','nfib','n_metal','ray_sqrt','scatter_num',...
    'dx_num','I_min_var','sma_num','max_scatter_ang','bounce_var','raysc_coeff','sma_fill_length','water_status'];
it_num=max(lengths_it);
index_it=find(lengths_it>1);%indexes of iterable variable(s)
legend_Main=" ";
leg_Main=" ";
%if/else statement for use in displaying data about the iteration/fiber
%bundle
if ~isempty(index_it)%if there are iterable variables%% TO DO generalize
    if length(index_it)>1%if there are multiple variables with multiple iterations make the set up
        % % yes_itname=2;
        % % it_num=sum((lengths_it(index_it)));
        % % Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
        % % x_len=ones(it_num,1).*SystemParam.xlen;%;reshape(Len,1,[]);
        % % scat_num=ones(it_num,1).*SystemParam.scatnum;%;reshape(Len,1,[]);
        % % sma_fill_length=ones(it_num,1).*SystemParam.SMA_fill;%;reshape(Len,1,[]);
        % % ang=ones(it_num,1).*SystemParam.scat_ang_max;%;reshape(Len,1,[]);
        % % rco=ones(it_num,1).*SystemParam.r_coeff;%;reshape(Len,1,[]);
        % % fibN=ones(it_num,1).*SystemParam.n1;%;reshape(Len,1,[]);
        % % metalN=ones(it_num,1).*SystemParam.n_metal;%;reshape(Len,1,[]);
        % % water_status=ones(it_num,1).*SystemParam.waterInterface;
        % % startind=1;
        % % if length(nfib)>1
        % %     endind=startind+length(nfib)-1;
        % %     fibN(startind:endind)=nfib;
        % %     nfib=fibN;
        % %     startind=endind+1;
        % % end
        % % if length(scatter_num)>1
        % %     endind=startind+length(scatter_num)-1;
        % %     scat_num(startind:endind)=scatter_num;
        % %     scatter_num=scat_num;
        % %     startind=endind+1;
        % % end
        % % if length(raysc_coeff)>1
        % %     endind=startind+length(raysc_coeff)-1;
        % %     rco(startind:endind)=raysc_coeff;
        % %     raysc_coeff=rco;
        % %     startind=endind+1;
        % % end
        % % if length(max_scatter_ang)>1
        % %     endind=startind+length(max_scatter_ang)-1;
        % %     ang(startind:endind)=max_scatter_ang;
        % %     max_scatter_ang=ang;
        % %     startind=endind+1;
        % % end
        % % if length(n_metal)>1
        % %     endind=startind+length(n_metal)-1;
        % %     metalN(startind:endind)=n_metal;
        % %     n_metal=metalN;
        % %     startind=endind+1;
        % % end
        % % 
        % % ang_denom=(max_scatter_ang./pi).^-1;
        % % leg_Main=cell(1,it_num);
        % % for i=1:it_num
        % %     leg_Main(1,i)={sprintf('water status; %.3g, scatter num; %.3g, maximum scatter angle denom; %.4g, rcoeff; %.3g, Imaginary RI; %.1e, Metal RI; %.3g, SMA fill amt; %.3g, fib length; %.3g',water_status(i),scatter_num(i),ang_denom(i),raysc_coeff(i),(imag(nfib(i))),real(n_metal(i)),sma_fill_length(i),x_len(i))};
        % % end
        % % legend_Main=string(leg_Main);
                yes_itname=2;
        it_num=prod(lengths_it(index_it));
        Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
        %[wat_stat,scat_num,ang,metalN,rco,fibN,sma_fill,Len]=ndgrid(water_status,scatter_num,max_scatter_ang,n_metal,raysc_coeff,nfib,sma_fill_length,x_len);
        [wat_stat,ang,rco,fibN,metalN,scat_num,Len,sma_fill]=ndgrid(water_status,max_scatter_ang,raysc_coeff,nfib,n_metal,scatter_num,x_len,sma_fill_length);

        %[ang,rco,fibN,metalN,scat_num]=ngrid(max_scatter_ang,raysc_coeff,nfib,n_metal,scatter_num);
        water_status=reshape(wat_stat,1,[]);
        x_len=reshape(Len,1,[]);
        scatter_num=reshape(scat_num,1,[]);
        max_scatter_ang=reshape(ang,1,[]);
        sma_fill_length=reshape(sma_fill,1,[]);
        ang_denom=(max_scatter_ang./pi).^-1;
        raysc_coeff=reshape(rco,1,[]);
        nfib=reshape(fibN,1,[]);
        n_metal=reshape(metalN,1,[]);
        leg_Main=cell(1,it_num);
        for i=1:it_num
            leg_Main(1,i)={sprintf('water status; %.3g, scatter num; %.3g, maximum scatter angle denom; %.4g, rcoeff; %.3g, Imaginary RI; %.1e, Metal RI; %.3g, SMA fill amt; %.3g, fib length; %.3g',water_status(i),scatter_num(i),ang_denom(i),raysc_coeff(i),(imag(nfib(i))),real(n_metal(i)),sma_fill_length(i),x_len(i))};
        end
        legend_Main=string(leg_Main);
%         error('More than one variable is being iterated, please fix and re-try')
    elseif SystemParam.num_fib>1 %if there is more than one fiber in the bundle and we're iterating (this should be a capability we have eventually though)
        yes_itname=0;
        error('You are trying to iterate different variables for a bundle of fibers, please fix and re-try')
    else%only one variable with multiple iteration. Determine which variable is being iterated
        yes_itname=1;%name of the iteration needs to be incorporated into data display
        switch index_it
            case 1
                it_name='Radius of Fiber (/mum)';
                Title_Main=sprintf('Irradiance along the Fiber Length, Varying Fiber Diameter \n %d Fiber(s), %d cm long, %d mm coupling distance',SystemParam.num_fib, xlen/10^4,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('%d (/mum) Diameter',rfiber(i)*2)};
                    end
                legend_Main=string(leg_Main);
            case 2
                it_name='Length of Fiber (/mum)';
                Title_main=sprintf('Irradiance along the Fiber Length, Varying Fiber Length \n %d Fiber(s),  %d [um] OD , %d mm coupling distance',SystemParam.num_fib,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('%d (cm] Length',x_len(i)/10^4)};
                    end
                legend_Main=string(leg_Main);
            case 3
                it_name='LED Distance from Fiber (mm)';
                Title_Main=sprintf('Irradiance along the Fiber Length, Varying Coupling Distance \n %d Fiber(s), %d cm long, %d [um] OD',SystemParam.num_fib, xlen/10^4,rfiber*2);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('%d (mm) Coupling Distance',led_dist(i))};
                    end
                legend_Main=string(leg_Main);
                %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 4
                it_name='imaginary part of complex refractive index';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('Imaginary RI: %d',imag(nfib(i)))};
                    end
                legend_Main=string(leg_Main);
                %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             case 5
                it_name=' metal refractive index';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('Imaginary metal RI: %d',n_metal(i))};
                    end
                legend_Main=string(leg_Main);
                %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             case 6
                it_name='Ray fidelity number';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('Ray fidelity number: %d',ray_sqrt(i))};
                    end
                legend_Main=string(leg_Main);
            case 7
                it_name='Scatter fidelity number';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('Scatter fidelity number: %d',scatter_num(i))};
                    end
                legend_Main=string(leg_Main);
            case 8
                it_name='dx fidelity';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('dx fidelity: %d [um]',dx_num(i))};
                    end
                legend_Main=string(leg_Main);
            case 9
                it_name='Minimum Intensity Tracked';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('Minimum Intensity Tracked: %d [uW/cm2]',I_min_var(i))};
                    end
                legend_Main=string(leg_Main);
            case 10
                it_name='number of housing bounces';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('housing bounce number: %d',sma_num(i))};
                    end
                legend_Main=string(leg_Main);
            case 11
                it_name='maximum scatter angle';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('maximum scatter angle: %d [rad]',max_scatter_ang(i))};
                    end
                legend_Main=string(leg_Main);
            case 12
                it_name='number of fiber bounces';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('fiber bounce number: %d',bounce_var(i))};
                    end
                legend_Main=string(leg_Main);
            case 13
                it_name='portion of light going to rayleigh scattering';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('rayleigh light %: %d',raysc_coeff(i))};
                    end
                legend_Main=string(leg_Main);
            case 14
                it_name='length to which the sma connector is filled with cytop';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('sma fill length%: %d [cm]',sma_fill_length(i))};
                    end
                legend_Main=string(leg_Main);
           case 15
                it_name='whether water is present or not';
                Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
                leg_Main=cell(1,lengths_it(index_it));
                    for i=1:lengths_it(index_it)
                        leg_Main(1,i)={sprintf('water medium (true or false)%: %d [cm]',water_status(i))};
                    end
                legend_Main=string(leg_Main);
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
        legend_Main=' ';
    end
    legend_Main=' ';
    Title_Main=sprintf('Irradiance along the Fiber Length \n %d Fiber(s), %d cm long, %d [um] OD , %d mm coupling distance',SystemParam.num_fib, xlen/10^4,rfiber*2,led_dist);
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


%%%%%%EACH ITERATION OF THE FIBER
for iteration=1:it_num
    r_fiber=SystemParam.r_fiber;
    %%%%%%%reset parameters to be reflective of the iteration number
    if yes_itname==1%if there are iterable items
        if index_it==1
            SystemParam.r_fiber=rfiber(iteration);
            r_fiber=SystemParam.r_fiber;
        elseif index_it==2
            SystemParam.xlen=xlen(iteration);
        elseif index_it==3
            SystemParam.led_distance=led_dist(iteration);
            %include cytop and/or NP variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif index_it==4
            SystemParam.n1=nfib(iteration);
        elseif index_it==5
            SystemParam.n_metal=n_metal(iteration);
        elseif index_it==6
            SystemParam.angnum=ray_sqrt(iteration);               %number of angles of led emission (must be less than raynum, below)
            SystemParam.raynum=SystemParam.angnum^2;  %total number of rays of led emission
            SystemParam.scatnum=ray_sqrt(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==7
            SystemParam.scatnum=scatter_num(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==8
            SystemParam.cont_dx=dx_num(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==9
            SystemParam.I_min=I_min_var(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==10
            SystemParam.housing_bounce=sma_num(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==11
            SystemParam.scat_ang_max=max_scatter_ang(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==12
            SystemParam.maxbounce=bounce_var(iteration);             %maximum number of rays produced by the scattering during end reflect
        elseif index_it==13
            r_coeff=raysc_coeff(iteration);
        elseif index_it==14
           SystemParam.SMA_fill= sma_fill_length(iteration);
        elseif index_it==15
            SystemParam.waterInterface=water_status(iteration);
        end
    elseif yes_itname==2 %to do generalize
        SystemParam.scat_ang_max=max_scatter_ang(iteration); 
        SystemParam.r_coeff=raysc_coeff(iteration);
        SystemParam.n1=nfib(iteration);
        SystemParam.n_metal=n_metal(iteration);
        SystemParam.waterInterface=water_status(iteration);
        SystemParam.scatnum=scatter_num(iteration);             %maximum number of rays produced by the scattering during end reflect
        SystemParam.SMA_fill=sma_fill_length(iteration);
        SystemParam.xlen=x_len(iteration);
        xlen=SystemParam.xlen;
        if x_len(iteration)>=48.5*10^4
            SystemParam.sealed_SMA=1;
        else
            SystemParam.sealed_SMA=0;
        end
    end

%calculate parameters that may need to use the changed iteration system
%param
    SystemParam.alpha=(10^-6)*(imag(SystemParam.n1)*4*pi/(SystemParam.uv_wavelength*10^-9))*10/log(10);%db/um overall absorption coefficient
% SystemParam.alpha_uv=alpha*(1-r_coeff);%db/um %ub absorption coeff
% SystemParam.alpha_r=alpha*r_coeff;%db/um %rayleigh scattering absorption coefficient

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
    LED_d=SystemParam.led_distance;%(um)%LED distance in um for Entering Light
    [Ray_X,Ray_Y,alpha_ang,beta_ang,Power_mat] = OutputLED3D(SystemParam);
    powercheck=sum(Power_mat,'all')*(SystemParam.LEDd*10^-1)^2;% make sure power isn't greater than the total LED power
    %determining light rays entering each fiber
    [circles,combocircles] = fiberBundle3D(SystemParam.num_fib,r_fiber*2);
    %[Ent_Int,Theta,y0,perc_hit,perc_ent,Entering_Int,Entering_angle,Entering_X,Entering_Y,incoming_int,incoming_ang] = enterLight3D_AINsubstrate(SystemParam,r_fiber, LED_d,alpha_ang,beta_ang,Ray_X,Ray_Y,Intensity_mat,circles);
    [Ent_Int,Theta,y0,perc_hit,perc_ent,Entering_Int,Entering_angle,Entering_X,Entering_Y,incoming_int,incoming_ang] = enterLight3D(SystemParam,r_fiber, LED_d,alpha_ang,beta_ang,Ray_X,Ray_Y,Power_mat,circles);
    %compressed to 2D
%    Ent_Int=Ent_Int;
    [a,b,c]=size(Ent_Int); %[num ang changes,num y locations, number fibers]

   
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
        aa_Rand(1).Pow_enter=zeros(1,aa_lim);
        aa_Rand(1).transmitted=zeros(1,aa_lim);
        aa_Rand(1).pow_side=zeros(1,aa_lim);
        aa_Rand(1).pow_side_use=zeros(1,aa_lim);
        aa_Rand(1).SMA_pow_side=zeros(1,aa_lim);
        aa_Rand(1).pow_side_waterst=zeros(1,aa_lim);
        aa_Rand(1).absorbed=zeros(1,aa_lim);
        aa_Rand(1).backscat=zeros(1,aa_lim);
        aa_Rand(1).SMAabs=zeros(1,aa_lim);
        aa_Rand(1).approxpow=zeros(1,aa_lim);
        aa_Rand(1).approxpow_dif=zeros(1,aa_lim);
        aa_Rand(1).approxpow_pos=zeros(1,aa_lim);
        aa_Rand(1).b2hpow=zeros(1,aa_lim);
        aa_Rand(1).cutoffpow=zeros(1,aa_lim);
        aa_Rand(1).remaininglosses=zeros(1,aa_lim);
        aa_Rand(1).UC=zeros(1,aa_lim);aa_Rand(2).UC=zeros(1,aa_lim);
        aa_Rand(1).RatioIsIt=zeros(1,aa_lim);aa_Rand(2).RatioIsIt=zeros(1,aa_lim);
        aa_Rand(1).Y=cell(aa_lim,1);aa_Rand(2).Y=cell(aa_lim,1);
%         aa_Rand(1).meas_inten=cell(1,aa_lim); %length of all the rays used
%         aa_Rand(1).meas_points=cell(1,aa_lim);%length of all the rays used
        aa_Rand(1).meas_inten=cell((a*b),aa_lim); %length of all the rays used
        aa_Rand(1).meas_points=cell((a*b),aa_lim);%length of all the rays used
        %empty storage vector setup for each individual ray 
        %possibly will need to include dimension for the scatter cone later
        %(aa_lim,a,b,gg_max)
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
                   if Ent_Int(xx,yy,h)==0%dont need to bother iterating through no intensity
                       continue
                   end
                    %measurement storage vectors
                    meas.points= zeros(2*10e6,2);
                    meas.inten = zeros(2*10e6,1);
                    meas.counter=1;
                    meas.sum=0;
                    
                    Global_Index=[iteration,h,aa,xx,yy];%record the index values, may need to include a gg if scatter cone at the front
                    Y_sum=zeros(size(Xvec));
                    %create a scatter cone
                    if SystemParam.scat_front==1
                        [I_scatter,Theta_enter]=scatter_cone(SystemParam,Theta(xx,yy,h),1);
                        I_enter=I_scatter.*Ent_Int(xx,yy,h);
                        gg_max=length(I_scatter);
                    else
                        I_enter=Ent_Int(xx,yy,h);
                        Theta_enter=Theta(xx,yy,h);
                        gg_max=1;
                    end
                        %tracking and summing values for each ray
                        %temporary storage of data
                        for gg=1:gg_max
                            %initialize the entering light values
                            I_ent=I_enter(gg);
                            Theta_ent=Theta_enter(gg);
                            %preallocate empty vectors with limits for the  while loop
                            max_loop=1+(SystemParam.maxbounce*SystemParam.scatnum);%max number of loops to run
                            
                            I_in=zeros(1,(max_loop));
                            I_intrack=zeros(1,(max_loop));
                            P_in=zeros((max_loop),2);
                            Theta_in=zeros(1,(max_loop));
                            Direction_in=zeros(1,(max_loop));
                            %assign starting while loop values
                            ray.Pow_enter(aa,xx,yy)=I_ent+ray.Pow_enter(aa,xx,yy);
                            I_in(1,1)=I_ent;
                            Theta_in(1,1)=Theta_ent;
                            P_in(1,1:2)=[0,y0(xx,yy,h)];
                            Direction_in(1,1)=1;
                            %assign initial index values
                            while_num=1;
                            bounce_num=1;
                            st_index=2;%starting index to allocate the first scattered ray we're tracking
                            pow_check=0;
                            MEAS0=sum(meas.inten);
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
                                %                     disp('Traveling')
                                Pow_out=[meas_dif,I_1,IT.backi,IT.absorbi,IT.housi,IT.cutoffi,IT.approxi];
                                
% %                                 [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_0,Pow_out,SystemParam,'Traveling',0,Global_Index,Tally);
% %                                 ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
% %                                 ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
                                
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
                                %                     disp('while loop step 0 to step 1')
                                Pow_Out=[meas_dif_1,I_1,b2h_1,cutoff_1,trans_1,abs_1,smaabs_1,approx_1,back_1];
                                %                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
% %                                 [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop step 0 to step 1',0,Global_Index,Tally);
                                
                                cutoffstart=ray.cutoffpow(aa,xx,yy);
                                meas_start2=sum(meas.inten);%initial total measured before entering into the end reflect function
                                %%%%%%%end reflect stuff%%%%%%%%%
                                %               disp('sum forward scatter')
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
                                            ray.cutoffpow(aa,xx,yy)=sum(I_scat(sort_scat(highest_index+1:end)))+ray.cutoffpow(aa,xx,yy);
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
                                cut_off=ray.cutoffpow(aa,xx,yy)-cutoffstart;
                                %                     disp('end reflect')
                                Pow_Out=[meas_dif,sum_forwardscat,cut_off,I_trans,b2h];
                                dif_endref=I_1-sum(Pow_Out);
% %                                 [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_1,Pow_Out,SystemParam,'End Reflect',0,Global_Index,Tally);
% %                                 ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
% %                                 ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
                                
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
                                %                     disp('while loop')
                                %                     start=[I_0,meas_start0,b2h0,cutoff0,trans0,abs0,smaabs0,approx0]
                                Pow_out=[sum_forwardscat,meas_dif_final,b2h_final,cutoff_final,trans_final,abs_final,smaabs_final,approx_final,back_final];
                                %                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
% %                                 [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop',0,Global_Index,Tally);
                                
                            end
                            %summing all of the measured side emitted power
                        end
                ray.pow_side(aa,xx,yy)=sum(meas.inten)-MEAS0;
%                 disp('total while loop')
                total_pow_use=[ray.pow_side(aa,xx,yy),ray.transmitted(aa,xx,yy),ray.SMAabs(aa,xx,yy),ray.b2hpow(aa,xx,yy),ray.backscat(aa,xx,yy),ray.absorbed(aa,xx,yy),ray.cutoffpow(aa,xx,yy)];
                
                %tracking differences
% %                 [Tally,Differenceamount,Diffamountpos] =  DifTrack(ray.Pow_enter(aa,xx,yy),total_pow_use,SystemParam,'total while loop',0,Global_Index,Tally);
% %                 ray.remaininglosses(aa,xx,yy)=Differenceamount;
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
                       [XVEC,Y, lengthPlot,pow_side_total] = Bins032725(meas_point_aa, meas_inten_aa, SystemParam.division, xlen, r_fiber,SystemParam);
                    use_index=find(meas_point_aa(:,1)>=(2.5*10^4));%indexes of all of the measurement points that correspond to measurable light after the sma
                    water_st_index=find(meas_point_aa(:,1)>=SystemParam.waterstart);
                    SMA_index=find(meas_point_aa(:,1)<=SystemParam.SMA_totallength);
                    pow_side_use=sum(meas_inten_aa(use_index));%*(2*pi*(SystemParam.r_fib*10^-4)*((SystemParam.;
                    SMA_pow_side=sum(meas_inten_aa(SMA_index));%*;
                    pow_side_waterst=sum(meas_inten_aa(water_st_index));
                    dAwater=((SystemParam.xlen*10^-4)-(SystemParam.waterstart*10^-4))*2*pi*(SystemParam.r_fiber*10^-4);
                    record_pow=pow_side_waterst*dAwater
                    %performance metric parameters
                    aa_Rand(1).Y(aa,1)=Y;
                    Ylocal=cell2mat(Y);%local version of the I(x) vector
                    aa_Rand(1).UC(1,aa)=Ylocal(i6)/Ylocal(i1);
                else
                    incr=SystemParam.division*10^-4;
                    if SystemParam.SMA==1
                    num_max = floor(((xlen/10e3)-(SystemParam.SMA_totallength/10e3))/(incr))+1;%dividing length of fiber by the increments, then adding one to have data at each end of bin
                
                    XVEC = zeros(1,num_max+2);    
                %first two measurement points will be within the SMA flush length
                Inc1=SystemParam.SMA_flushlength*10^-4;
                Inc2=(SystemParam.SMA_totallength-SystemParam.SMA_flushlength)*10^-4;
                
                %first two points have irregular indexes if a part of the sma connector
                XVEC(1)=-(Inc2+Inc1);
                XVEC(2)=-Inc1;
                XVEC(3:end)=0:incr:(incr*(num_max-1));
                    else
                      num_max = floor((xlen/10e3)/(incr))+1;%dividing length of fiber by the increments, then adding one to have data at each end of bin

                        XVEC=0:incr:(incr*(num_max-1));
                    end
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
        trans_pow=aa_Rand(1).transmitted(1,aa).*pi*(SystemParam.r_fiber*10^-4)^2

        aa_Rand(1).RatioIsIt(aa,xx,yy)=aa_Rand(1).pow_side(1,aa)/aa_Rand(1).transmitted(1,aa);
        aa_Rand(1).pow_side(1,aa)=pow_side_total;
        aa_Rand(1).pow_side_use(1,aa)=pow_side_use;
        aa_Rand(1).SMA_pow_side(1,aa)=SMA_pow_side;
        aa_Rand(1).pow_side_waterst(1,aa)=pow_side_waterst;
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
        FibIt(1).pow_side(iteration,h)=mean(aa_Rand(1).pow_side(1,:));
        FibIt(1).pow_side_use(iteration,h)=mean(aa_Rand(1).pow_side_use(1,:));
        FibIt(1).SMA_pow_side(iteration,h)=mean(aa_Rand(1).SMA_pow_side(1,:));
        FibIt(1).pow_side_waterst(iteration,h)=mean(aa_Rand(1).pow_side_waterst(1,:));
        FibIt(1).absorbed(iteration,h)=mean(aa_Rand(1).absorbed(1,:));
        FibIt(1).backscat(iteration,h)=mean(aa_Rand(1).backscat(1,:));
        FibIt(1).SMAabs(iteration,h)=mean(aa_Rand(1).SMAabs(1,:));
        FibIt(1).approxpow(iteration,h)=mean(aa_Rand(1).approxpow(1,:));
        FibIt(1).approxpow_dif(iteration,h)=mean(aa_Rand(1).approxpow_dif(1,:));
        FibIt(1).approxpow_pos(iteration,h)=mean(aa_Rand(1).approxpow_pos(1,:));
        FibIt(1).b2hpow(iteration,h)=mean(aa_Rand(1).b2hpow(1,:));
        FibIt(1).cutoffpow(iteration,h)=mean(aa_Rand(1).cutoffpow(1,:));
        FibIt(1).remaininglosses(iteration,h)=mean(aa_Rand(1).remaininglosses);

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
        FibIt(2).remaininglosses(iteration,h)=std(aa_Rand(1).remaininglosses);

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
                

%          %storing the larger data for maybe checking the trends out with
%          SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
%          SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
%          SumVal.transmitted(iteration,h)={ray.transmitted};
%          SumVal.pow_side(iteration,h)={ray.pow_side};
%          SumVal.pow_side_use(iteration,h)={ray.pow_side_use};
%          SumVal.SMA_pow_side(iteration,h)={ray.SMA_pow_side};
%          SumVal.pow_side_waterst(iteration,h)={ray.pow_side_waterst};
%          SumVal.absorbed(iteration,h)={ray.absorbed};
%          SumVal.backscat(iteration,h)={ray.backscat};
%          SumVal.SMAabs(iteration,h)={ray.SMAabs};
%          SumVal.approxpow_dif(iteration,h)={ray.approxpow_dif};%sum of power differences due to difference in approximations
%          SumVal.approxpow_pos(iteration,h)={ray.approxpow_pos}; %absolute value of all of the power difference approximationsSumVal.b2hpow(iteration,h)={ray.;
%          SumVal.cutoffpow(iteration,h)={ray.cutoffpow};
%          SumVal.remaininglosses(iteration,h)={aa_Rand(1).remaininglosses};
%          SumVal.UC(iteration,h)={ray.UC};
%          SumVal.RatioIsIt(iteration,h)={ray.RatioIsIt};
%          SumVal.Y(iteration,h)={ray.Y};
         
         %plotting the pie chart tracking stuff
         %percent of each thing for power balance
figure(1)
Power_Track=[aa_Rand(1).pow_side(1,aa),aa_Rand(1).transmitted(1,aa),aa_Rand(1).SMAabs(1,aa),aa_Rand(1).b2hpow(1,aa),aa_Rand(1).backscat(1,aa),aa_Rand(1).absorbed(1,aa),aa_Rand(1).cutoffpow(1,1),aa_Rand(1).remaininglosses(1,1)];
Power_Track_lab=["pow_side" "transmitted" "SMAabs" "b2hpow" "backscat" "absorbed" "cutoffpow" "remaininglosses"];
pie(Power_Track,Power_Track_lab)
% 
% figure(2)
% Power_Track_comp=[aa_Rand(1).pow_side_use(1,aa),aa_Rand(1).SMA_pow_side(1,aa),aa_Rand(1).transmitted(1,1),aa_Rand(1).SMAabs(1,1),aa_Rand(1).b2hpow(1,1),aa_Rand(1).backscat(1,1),aa_Rand(1).absorbed(1,1),aa_Rand(1).cutoffpow(1,1),aa_Rand(1).remaininglosses(1,1)];
% Power_Track_comp_lab=["pow_side_use" "SMA_pow_side" "transmitted" "SMAabs" "b2hpow" "backscat" "absorbed" "cutoffpow" "remaininglosses"];
% pie(Power_Track_comp,Power_Track_comp_lab)
% 
% %plotting the I(x) distribution data for each iteration
% figure(3)
% %title and legend defined earlier lines 95 is through 140ish
% Y_current=cell2mat(FibIt(1).Y(iteration,h));%current I(x) within this iteration and fiber
% plot(XVEC,Y_current)
% title(Title_Main) %defined earlier 
% hold on
% 
% 

% % 
% % %%%%% this is all just extracting the data from the Tally struct to find where the biggest differences in coming from for trouble shooting
% % %%%% i didn't comment this well nor is this probably the best method to
% % %%%% extract the data. but it works for now
% % %figuring out where the gaps are, sort the strings of the difference values
% % [m1,m2,m3,m4]=size((Tally.big_dif_fun_main));%[iteration, fiber number, aa limit, number of times thing got assigned to the tally]
% % [f1,f2,f3,f4,f5,f6]=size((Tally.big_dif_fun_func));%[iteration, fiber number, aa limit,xx,yy, number of times thing got assigned to the tally]
% % 
% % maincount=numel(Tally.big_dif_count_meas_main);%total number of times we count a "big difference" in the main code
% % funcount=numel(Tally.big_dif_count_meas_func);%total number of times we count a "big difference" within each function
% % IO_maincount=numel(Tally.IO_count_meas_main);
% % IO_funcount=numel(Tally.IO_count_meas_func);
% % 
% % bigmain=strings(1,maincount);%empty vector for identifying the part of the code that has the big difference
% % dif_main=zeros(1,maincount);%empty vector for recording the big difference within the main
% % dif_fun=zeros(1,funcount);%empty vector for identifying the function that has the big difference
% % bigfun=strings(1,funcount);%empty vector for recording the big difference within the function
% % cmain=1;%index
% % cfun=1;%index for the main code
% % %for the in and out check
% % IO_bigmain=strings(1,IO_maincount);%empty vector for identifying the part of the code that has the big difference
% % IO_dif_main=zeros(1,IO_maincount);%empty vector for recording the big difference within the main
% % IO_dif_fun=zeros(1,IO_funcount);%empty vector for identifying the function that has the big difference
% % IO_bigfun=strings(1,IO_funcount);%empty vector for recording the big difference within the function
% % IO_cmain=1;%index
% % IO_cfun=1;%index for the main code
% % 
% % for i=1:m1 %iteration number
% %     for j=1:m2 %fiber number
% %         for k=1:m3 %aa number
% %             %set up emtry string vector
% %            count_inst=Tally.big_dif_count_main(i,j,k);%find maximum index where there is a location in main that had
% %            IO_count_inst=Tally.IO_count_main(i,j,k);%find maximum index where there is a location in main that had
% % 
% %            %a big difference
% %             for l=1:(count_inst)
% %                 if ischar(Tally.big_dif_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
% %                    
% %                     bigmain(1,cmain)=Tally.big_dif_fun_main{i,j,k,l}{1,1};%assign the string value
% %                     dif_main(1,cmain)=cell2mat(Tally.big_dif_count_meas_main{i,j,k,l}); %assign the difference value
% %                     cmain=cmain+1; %update the index
% %                 end
% %             end
% %             for l=1:(IO_count_inst)
% %                 if ischar(Tally.IO_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
% %                    
% %                     IO_bigmain(1,IO_cmain)=Tally.IO_fun_main{i,j,k,l}{1,1};%assign the string value
% %                     IO_dif_main(1,IO_cmain)=cell2mat(Tally.IO_count_meas_main{i,j,k,l}); %assign the difference value
% %                     IO_cmain=IO_cmain+1; %update the index
% %                 end
% %             end
% %             %assign
% %             for l=1:f4
% %                 for m=1:f5
% %                     %find each location where there is a function that had
% %                     %a big difference
% %                     count_inst=find(Tally.big_dif_fun_func{i,j,k,l,m});%
% %                     IO_count_inst=Tally.IO_count_func(i,j,k,l,m);
% % 
% %                     for o=1:length(count_inst)
% %                         n=count_inst(o);%the index
% %                             %double check the difference function isn't
% %                             %empty
% %                             if isempty(Tally.big_dif_fun_func{i,j,k,l,m,n})
% %                                 fun_str_inst="";
% %                             else
% %                                 fun_str_inst=Tally.big_dif_fun_func{i,j,k,l,m,n};
% %                                  dif_fun(1,cfun)=Tally.big_dif_count_meas_func{i,j,k,l,m,n};
% %                                  bigfun(1,cfun)=fun_str_inst;
% %                                  cfun=cfun+1;
% % 
% %                             end
% %                     end
% %                     if IO_count_inst~=0
% %                         for o=1:length(IO_count_inst)
% %                             n=IO_count_inst(o);%the index
% %                                 %double check the difference function isn't
% %                                 %empty
% %                                 if isempty(Tally.IO_fun_func{i,j,k,l,m,n})
% %                                     fun_str_inst="";
% %                                 else
% %                                     fun_str_inst=Tally.IO_fun_func{i,j,k,l,m,n};
% %                                     IO_dif_fun(1,IO_cfun)=Tally.IO_count_meas_func{i,j,k,l,m,n};
% %                                     IO_bigfun(1,IO_cfun)=fun_str_inst;
% %                                     IO_cfun=IO_cfun+1;
% % 
% %                                 end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % idfun=find(bigfun~="");
% %     if sum(dif_fun,'all')==sum(dif_fun(idfun),'all')
% %         bigfun=bigfun(idfun);
% %         fun_cats=unique(bigfun);%create a category of each unique function reported with a big dif
% %         trouble_functions=categorical(bigfun,fun_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % idmain=find(bigmain~="");
% %     if sum(dif_main,'all')==sum(dif_main(idmain),'all')
% %         bigmain=bigmain(idmain);
% %         main_cats=unique(bigmain);%create a category of each unique function reported with a big dif
% %         trouble_main=categorical(bigmain,main_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % %main_cats=unique(bigmain);%create a category of each unique function reported with a big dif
% % summary(trouble_functions)%get a summary of all of the trouble functions
% % %trouble_main=categorical(bigmain,main_cats);%main_categorical);
% % summary(trouble_main) %get a summary of all of the trouble areas in the main code
% % 
% % IO_idfun=find(IO_bigfun~="");
% %     if sum(IO_dif_fun,'all')==sum(IO_dif_fun(IO_idfun),'all')
% %         IO_bigfun=IO_bigfun(IO_idfun);
% %         IO_fun_cats=unique(IO_bigfun);%create a category of each unique function reported with a big dif
% %         IO_trouble_functions=categorical(IO_bigfun,IO_fun_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % IO_idmain=find(IO_bigmain~="");%    ;
% % IO_main_cats=unique(IO_bigmain(IO_idmain));%create a category of each unique function reported with a big dif
% % summary(IO_trouble_functions);%get a summary of all of the trouble functions
% % IO_trouble_main=categorical(IO_bigmain,IO_main_cats);%main_categorical);
% % summary(IO_trouble_main) %get a summary of all of the trouble areas in the main code
% % 
% % %create a table with summary values of  issues caused by each
% % %function
% % %empty cell for all of the categories of each big_dif item
% % dif_f=cell(length(fun_cats),1);
% % dif_m=cell(length(main_cats),1);
% % IO_dif_f=cell(length(fun_cats),1);
% % IO_dif_m=cell(length(main_cats),1);
% % %empty vectory for all of the categories of each big_dif item
% % dif_contrib_m=zeros(length(main_cats),2);
% % dif_contrib_f=zeros(length(fun_cats),2);
% % 
% % IO_dif_contrib_m=zeros(length(IO_main_cats),2);
% % IO_dif_contrib_f=zeros(length(IO_fun_cats),2);
% % for i=1:length(main_cats)%for each categorical in the main function
% %     dif_m(1,i)={dif_main(1,(trouble_main==main_cats(i)))};% max a cell 
% %     dif_contrib_m(i,1)=sum(cell2mat(dif_m(1,i)),'all');
% %     dif_contrib_m(i,2)=sum(abs(cell2mat(dif_m(1,i))),'all');
% % end
% % for i=1:length(fun_cats)
% %     dif_f(1,i)={dif_fun(1,(trouble_functions==fun_cats(i)))};
% %     dif_contrib_f(i,1)=sum(cell2mat(dif_f(1,i)),'all');
% %     dif_contrib_f(i,2)=sum(abs(cell2mat(dif_f(1,i))),'all');
% % end
% % 
% % for i=1:length(IO_main_cats)%for each categorical in the main function
% %     IO_dif_m(1,i)={IO_dif_main(1,(IO_trouble_main==IO_main_cats(i)))};% max a cell 
% %     IO_dif_contrib_m(i,1)=sum(cell2mat(IO_dif_m(1,i)),'all');
% %     IO_dif_contrib_m(i,2)=sum(abs(cell2mat(IO_dif_m(1,i))),'all');
% % end
% % for i=1:length(IO_fun_cats)
% %     IO_dif_f(1,i)={IO_dif_fun(1,(IO_trouble_functions==IO_fun_cats(i)))};
% %     IO_dif_contrib_f(i,1)=sum(cell2mat(IO_dif_f(1,i)),'all');
% %     IO_dif_contrib_f(i,2)=sum(abs(cell2mat(IO_dif_f(1,i))),'all');
% % end
% % 
% % 
% % Changed=dif_contrib_m(:,1);%difference value vector for the table
% % Total=dif_contrib_m(:,2);%absolute value vector for the table
% % Function=transpose(main_cats); %function
% % TabM=table(Changed,Total,'RowNames',Function)
% % Changed=dif_contrib_f(:,1);%difference value vector for the table
% % Total=dif_contrib_f(:,2);%absolute value vector for the table
% % Function=transpose(fun_cats);
% % TabF=table(Changed,Total,'RowNames',Function)
% % 
% % 
% % IO_Changed=IO_dif_contrib_m(:,1);%difference value vector for the table
% % IO_Total=IO_dif_contrib_m(:,2);%absolute value vector for the table
% % IO_Function=transpose(IO_main_cats); %function
% % IO_TabM=table(IO_Changed,IO_Total,'RowNames',IO_Function)
% % IO_Changed=IO_dif_contrib_f(:,1);%difference value vector for the table
% % IO_Total=IO_dif_contrib_f(:,2);%absolute value vector for the table
% % IO_Function=transpose(IO_fun_cats);
% % IO_TabF=table(IO_Changed,IO_Total,'RowNames',IO_Function)
% % 

%%%%%%writing the data to an excel file
    %%%%SAVING EACH ITERATION, each fiber DATA
    %saving
        if (writeToFile == 1)
            writeResults(SystemParam, description, filename, sheet_num, iteration, h, FibIt, c, legend_Main, aa_lim, LED_d, XVEC)
        end
    end%%%%%%%%%%
%  load gong.mat%gong to signal the code is done handel.mat %chorus sound to signal the code is done %
%  sound(y)

end
% if yes_itname==1%if we need a legend bc there's multiple iterations
%     legend(legend_Main)
% end
hold off
% load gong.mat%gong to signal the code is done handel.mat %chorus sound to signal the code is done %
% sound(y)

