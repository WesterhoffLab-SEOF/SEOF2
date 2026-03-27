function[I_return,P_return,V_return,IT,Tally,meas]=meas_SMA_check(I_in,I_transmit,I_reflect,Vtransmit,Vref,P,SystemParam,IT,Tally,meas,Global_Index)
%Function checks if the SMA connector is in the way before taking a
%measurement. if in the sma it runs a while loop that bounces a ray back into the fiber
%TO DO: include coating conditions for surface roughness model
intensityMin=SystemParam.intensityMin;
%standard outputs
P_return=P;
I_return=I_reflect;
V_return=Vref;

%if the transmitted light is functionally 0 (ie emitted into the sma
%connector flush part, or just very tiny to begin with)
if I_transmit<=intensityMin
    %skip any measurements
    %exit the code
    return
end

%check there is an SMA. can probably skip this
if SystemParam.SMA==0%if there isn't an sma connector then all transmitted light is measurable
    %measurements
    meas.inten(meas.counter)=I_transmit;%measure the transmitting light
    meas.points(meas.counter,:)=P;%measure the location
    meas.counter=meas.counter+1;%increase the counter position
    meas.sum=I_transmit+meas.sum;
    if meas.inten(meas.counter-1)==0
        error('measuring 0')
    end
    return
end
%if there is an sma connector, need to  set up for the sma connector
SMAFL=SystemParam.smaFlushLength;%=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SMA_L=SystemParam.smaTotalLength;%=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SMA_d=SystemParam.smaDiameter;%diameter of SMA connector [mm to um]
bounce_num=SystemParam.housingBounce;%maximum number of times I'm willing to let the ray bounce between the SMA and the fiber
iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);%indexes
%for things inside of the SMA connector
%parametric travel
dySMA=sign(Vtransmit(2))*((SMA_d/2)-SystemParam.fiberRadius);%distance in the direction of travel to what would be the sma edhge
tSMA=dySMA/Vtransmit(2);
dxSMA=tSMA*Vtransmit(1);
XSMA=P(1)+dxSMA;%x location at possible sma edge
%bounce =1;
%first the vtransmit isn't horizontal, or vertical
ratYtoX=abs(Vtransmit(2)/Vtransmit(1));
ratXtoY=abs(Vtransmit(1)/Vtransmit(2));

if ratYtoX==0 %if the transmitting ray is  entirely horizontl, not worried about bounce
    if P(1)>SMAFL %all reflecteed light is returned, but the transmitted light is considered either transmitted or returned to housing
        meas.inten(meas.counter)=I_transmit;%measure the transmitting light
        meas.points(meas.counter,:)=P;%measure the location
        meas.counter=meas.counter+1;%increase the counter position
        meas.sum=I_transmit+meas.sum;
        if meas.inten(meas.counter-1)==0
            display(I_transmit)
            display(meas.inten(meas.counter-1))
            error('measuring 0')
        end

    end
    return
elseif ratXtoY<=0 %if the transmitting ray is only vertical
    Tally.photon_vert_count(iteration,h,aa)=Tally.photon_vert_count(iteration,h,aa)+1;%counter goes up
    P_return=P;

    if P(1)<=SMAFL && P(1)>=0 %if it's within the flush length, lose the light

    else
        %record the measureable values. assume all light is side
        %emitted
        meas.points(meas.counter,:)=P;
        meas.inten(meas.counter,:)=I_transmit+I_reflect;
        meas.counter=meas.counter+1;
        meas.sum=I_transmit+I_reflect+meas.sum;
        if meas.inten(meas.counter-1)==0
            display(I_transmit)
            display(meas.inten(meas.counter-1))
            error('measuring 0')
        end
        I_return=I_reflect;%I_reflect;%have no light returning eventually

    end
    return
elseif P(1)>SMAFL && XSMA>SMA_L%if the impact point isn't within the SMA connector & the transmitting ray wont bounce off the SMA connector
    %can take a normal measurement and go back to the traveling function
    %measurements
    meas.inten(meas.counter)=I_transmit;%measure the transmitting light
    meas.points(meas.counter,:)=P;%measure the location
    meas.counter=meas.counter+1;%increase the counter position
    meas.sum=I_transmit+meas.sum;
    if meas.inten(meas.counter-1)==0
        display(I_transmit)
        display(meas.inten(meas.counter-1))
        error('measuring 0')
    end
    %set up for returning to the traveling function
    I_return=I_reflect;

    return
elseif P(1)>=0 && P(1)<=SMAFL %if the ray impacts within the flush part of the sma connector
    %everything is reflected, nothing is measured
    %             [ni,k,Tco]=medium_check(SystemParam,P(1));%quick function to figure what medium we're in + the relevant proerties
    I_return=I_reflect;
    return
elseif XSMA<SMAFL %if the light bounced outside of the flush part, but would hit a complicated part of the interior SMA geometry
    %assume it to be not bounce-able/unimportant to follow
    %take a normal measurement and reset
    %measurements
    meas.inten(meas.counter)=I_transmit;%measure the transmitting light
    meas.points(meas.counter,:)=P;%measure the location
    meas.counter=meas.counter+1;%increase the counter position
    meas.sum=I_transmit+meas.sum;
    if meas.inten(meas.counter-1)==0
        display(I_transmit)
        display(meas.inten(meas.counter-1))
        error('measuring 0')
    end
    %set up for returning to the traveling function
    P_return=P;
    I_return=I_reflect;
    return
else%if the light bounces in the sma (outside of the flush part), and the ray would hit the bounceable part of the sma connector
    %set up conditions and initial inputs for a while loop to track the
    %bouncing and losses
    dirVt=sign(Vtransmit(1));%calculate the direction of x travel
    Vbounce=Vtransmit;%this should stay the same the whole time (assumign flat surface on top)
    Pbounce=P;
    xbounce=P(1);%this
    bounce=1;
    I_bounce=I_transmit;
    I_bounceloss=0;%tally of how much loss there is
    nhat_metal=[0,-1*sign(Vbounce(2))];
    horz_surf=0;
    switch dirVt
        case -1%
            conditionSMA=(XSMA>SMAFL);
            condition=(xbounce>0);
        case 1
            conditionSMA=(XSMA<SMA_L);
            condition=(xbounce<SystemParam.xlen);
    end
    I_extra_return=zeros(1,bounce_num);
    P_extra_return=zeros(bounce_num,2);
    V_extra_return=zeros(bounce_num,2);


    while bounce<=bounce_num && condition && conditionSMA &&I_bounce>intensityMin
        %take measurement at the top using the I_bounce
        meas0=sum(meas.inten);%house0=IT.housi;
        meas.inten(meas.counter)=I_bounce;
        meas.points(meas.counter,:)=Pbounce;
        meas.counter=meas.counter+1;
        meas.sum=I_bounce+meas.sum;

        if meas.inten(meas.counter-1)==0
            display(I_transmit)
            display(meas.inten(meas.counter-1))
            error('measuring 0')
        end
        drSMA=norm(dxSMA,dySMA);%distance traveled to get to SMA surface
        xbounce1=(dxSMA)+Pbounce(1);%new xbounce
        %check what medium the ray travels through
        [ni,k,Tco]=medium_check(SystemParam,xbounce);%quick function to figure what medium we're in + the relevant proerties
        n_metal=SystemParam.n_metal;%
        %calc attenuation coeff over the travel dist to hit SMA
        kloss1=Tco*exp(-k*(drSMA*10^-4));
        I_bouncei1=kloss1*I_bounce;
        %rtrack atten loss
        I_bounceloss=I_bounce*(1-kloss1)+I_bounceloss;%loss amount
        %do reflection calc off of the metal
        [theta_im,theta_tm,theta_cm,~,~,~,V_reflectm,~] = Snells(Vbounce,nhat_metal,ni,n_metal,horz_surf,dirVt);
        %[I_reflect_m,I_lost_m]=FresnelEq(I_bouncei1,SystemParam,theta_im,theta_tm,theta_cm,theta_ihm,ni,n_metal,horz_surf);
        [I_reflect_m,I_lost_m]=FresnelEqSEOFv2(I_bouncei1,theta_im,theta_tm,theta_cm,ni,n_metal);

        %track amt lost to metal
        I_bounceloss=I_lost_m+I_bounceloss;%loss amount

        %rest the direction calculations for the travel distance back to the fiber
        %using reflected metal vel
        dySMA=sign(V_reflectm(2))*((SMA_d/2)-SystemParam.r_fiber);%distance in the direction of travel to what would be the sma edhge
        tSMA=dySMA/V_reflectm(2);
        dxSMA=tSMA*V_reflectm(1);
        drSMA=norm(dxSMA,dySMA);
        %reset jxbounce
        xbounce=xbounce1+dxSMA;
        Vreturn=V_reflectm;

        %determine the refractive indices && attenuation coefficients
        nt=SystemParam.n1;

        kloss=exp(-k*(drSMA*10^-4));


        I_bouncei=kloss*I_reflect_m;%input intensity into fresnel eq()
        %track atten loss
        I_bounceloss=I_reflect_m*(1-kloss)+I_bounceloss;%loss amount
        %define the surface
        nhat=[0,-1];
        dirVt=sign(Vreturn(1));
        horz_surf=1;
        %snells law for angles and velocity
        %fresnell equation
        [theta_i,theta_t,theta_c,~,~,V_transmit,V_reflect,~] = Snells(Vreturn,nhat,ni,nt,horz_surf,dirVt);

        %[I_reflect_sma,I_transmit_fib]=FresnelEq(I_bouncei,SystemParam,theta_i,theta_t,theta_c,theta_ih,ni,nt,horz_surf);
        [I_reflect_sma,I_transmit_fib]=FresnelEqSEOFv2(I_bouncei,theta_i,theta_t,theta_c,ni,nt);
        %check if there's a lot of differences
        %record the extra intensity coming back in
        I_extra_return(bounce)=I_transmit_fib;

        I_bounce=I_reflect_sma;
        Pbounce=[xbounce,P(2)];
        P_extra_return(bounce,:)=Pbounce;
        V_extra_return(bounce,:)=V_transmit;
        Vbounce=V_reflect;

        if I_bounce<SystemParam.intensityMin
            Tally.minPhotons_count(iteration,h,aa)=Tally.minPhotons_count(iteration,h,aa)+1;%counting all the times the model exits a function bc the direction is vertical
            break
        end
        %update for next while loop
        bounce=bounce+1;
        %

        %check the conditions
        switch dirVt
            case -1%
                conditionSMA=(XSMA>SMAFL);
                condition=(xbounce>0);
            case 1
                conditionSMA=(XSMA<SMA_L);
                condition=(xbounce<SystemParam.xlen);
        end
    end
    IT.cutoffi=IT.cutoffi+I_bounce;%cut off the last bit  of I_bounce
    IT.housi=I_bounceloss+IT.housi;
    I_return=I_reflect+sum(I_extra_return);

    %ASSUME ALL the returning light is concentrated at the same point with a
    %weighted average

    P_return(1)=((P(1).*I_reflect)+(dot(P_extra_return(:,1),I_extra_return)))./I_return;
    P_return(2)=P(2);
    V_return(1)=((Vref(1).*I_reflect)+(dot(V_extra_return(:,1),I_extra_return)))./I_return;
    V_return(2)=((Vref(2).*I_reflect)+(dot(V_extra_return(:,2),I_extra_return)))./I_return;

    return
end
end
