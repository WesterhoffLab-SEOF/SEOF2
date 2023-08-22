function[I_return,P_return,IT,Tally,meas]=meas_SMA_check(I_in,I_transmit,I_reflect,Vtransmit,P,SystemParam,IT,Tally,meas,Global_Index)
%Function checks if the SMA connector is in the way before taking a
%measurement
%measurement, then runs a while loop that bounces a ray back into the fiber
%include (update this to figure out how long to track a bouncing ray in a
%coating
%if it's in the right place
%%    Meas0=sum(meas.inten);House0=IT.housi;Absorb0=IT.absorbi;Back0=IT.backi;Cutoff0=IT.cutoffi;
I_min=SystemParam.I_min;
%check there is an SMA
if SystemParam.SMA==0%if there isn't an sma connector then all transmitted light is measurable
    %measurements
    meas.inten(meas.counter)=I_transmit;%measure the transmitting light
    meas.points(meas.counter,:)=P;%measure the location
    meas.counter=meas.counter+1;%increase the counter position
    meas.sum=I_transmit+meas.sum;
    if meas.inten(meas.counter-1)==0
        error('measuring 0')
    end
    %set up for returning to the traveling function
    P_return=P;
    I_return=I_reflect;
%%     [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten))],SystemParam,'No SMA w/in meas_SMA_check',1,Global_Index,Tally);
%%     [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'No SMA w/in meas_SMA_check',1,Global_Index,Tally);

    return
end
%then set up for the sma connector
metal_abs=SystemParam.metal_abs;%percent of light absorbed by hitting the SMA metal
SystemParam.SMA=1;%is there a SMA connector used to couple thhe fiber to the LED? 1 for yes, 0 for no
SMAFL=SystemParam.SMA_flushlength;%=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SMA_L=SystemParam.SMA_totallength;%=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SMA_d=SystemParam.SMA_diam;%diameter of SMA connector [mm to um]
bounce_num=SystemParam.housing_bounce;%maximum number of times I'm willing to let the ray bounce between the SMA and the fiber
iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);%indexes
%for things inside of the SMA connector
%parametric travel
dySMA=sign(Vtransmit(2))*((SMA_d/2)-SystemParam.r_fiber);%distance in the direction of travel to what would be the sma edhge
tSMA=dySMA/Vtransmit(2);
dxSMA=tSMA*Vtransmit(2);
XSMA=P(1)+dxSMA;%x location at possible sma edge
%bounce =1;
%first the vtransmit isn't horizontal, or vertical 
ratYtoX=abs(Vtransmit(2)/Vtransmit(1));
ratXtoY=abs(Vtransmit(1)/Vtransmit(2));

if ratYtoX==0 %if the transmitting ray is  entirely horizontl, not worried about bounce
    %check it's not within the flush part of the SMA connector
%%            b2h0=IT.b2hi;
%%        trans0=IT.transi;
    if P(1)>SMAFL %all reflecteed light is returned, but the transmitted light is considered either transmitted or returned to housing

        I_return=I_reflect;

        meas.inten(meas.counter)=I_transmit;%measure the transmitting light
        meas.points(meas.counter,:)=P;%measure the location
        meas.counter=meas.counter+1;%increase the counter position
        meas.sum=I_transmit+meas.sum;
        if meas.inten(meas.counter-1)==0
            error('measuring 0')
        end
    else %if it is within the sma connector, remove some of the absorbing light, can't measure
        I_return=I_reflect*(1-metal_abs);
        IT.housi=IT.housi+(I_reflect*metal_abs)+I_transmit;
    end
    P_return=P;
    %%disp(I_transmit)
    %%disp(I_reflect)

%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(sum(meas.inten)-Meas0),(IT.housi-House0),(IT.transi-trans0),(IT.b2hi-b2h0)],SystemParam,'horz ray w/in meas_SMA_check',1,Global_Index,Tally);
%%            [Tally,~,~]=inOutTrack(I_in,(sum(meas.inten)-Meas0),SystemParam,'horz ray w/in meas_SMA_check',1,Global_Index,Tally);

    return
elseif ratXtoY<=0 %if the transmitting ray is only vertical
        Tally.photon_vert_count(iteration,h,aa)=Tally.photon_vert_count(iteration,h,aa)+1;%counter goes up
    P_return=P;

        if P(1)<=SMAFL && P(1)>=0 %if it's within the flush length, lose the light
            %record the values
            %assume we lose all of the light into the housing, nothing is
            %measured
            IT.housi=I_transmit+I_reflect+IT.housi;
            I_return=0;


        else
            %record the measureable values. assume all light is side
            %emitted
            meas.points(meas.counter,:)=P;
            meas.inten(meas.counter,:)=I_transmit+I_reflect;
            meas.counter=meas.counter+1;
            meas.sum=I_transmit+I_reflectmeas.sum;
            if meas.inten(meas.counter-1)==0
                display(I_transmit)
                display(meas.inten(meas.counter-1))
                error('measuring 0')
            end
            I_return=I_reflect;%have no light returning eventually

        end
%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten)),(House0-IT.housi)],SystemParam,'vert ray w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'vert ray w/in meas_SMA_check',1,Global_Index,Tally);

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
    P_return=P;
    I_return=I_reflect;
%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten)),(House0-IT.housi)],SystemParam,'SMA1, normal, w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'SMA1, normal, w/in meas_SMA_check',1,Global_Index,Tally);

    return
elseif P(1)>=0 && P(1)<=SMAFL %if the ray impacts within the flush part of the sma connector
    %everything is reflected, nothing is measured
    I_return=I_reflect*(1-metal_abs);%mostly things reflect except for whats lost during
    IT.housi=(I_reflect*(metal_abs))+I_transmit+IT.housi;%amount absorbed by the metal
    P_return=P;
%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten)),(House0-IT.housi)],SystemParam,'SMA total loss w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'SMA total loss w/in meas_SMA_check',1,Global_Index,Tally);

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
        error('measuring 0')
    end
    %set up for returning to the traveling function
    P_return=P;
    I_return=I_reflect;
%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten)),(House0-IT.housi)],SystemParam,'SMA endbouncequick w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'SMA endbouncequick w/in meas_SMA_check',1,Global_Index,Tally);

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
    


    while bounce<=bounce_num && condition && conditionSMA &&I_bounce>I_min
            %take measurement at the top using the I_bounce
            meas0=sum(meas.inten);house0=IT.housi;
    meas.inten(meas.counter)=I_bounce;
    meas.points(meas.counter,:)=Pbounce;
    meas.counter=meas.counter+1;
    meas.sum=I_bounce+meas.sum;
%%    [Tally,~,~]=DifTrack(I_bounce,(meas0-sum(meas.inten)),SystemParam,'bounceloop measurement w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(meas0-sum(meas.inten)),SystemParam,'bounceloop measurement w/in meas_SMA_check',1,Global_Index,Tally);

    if meas.inten(meas.counter-1)==0
        display(I_transmit)
        display(meas.inten(meas.counter-1))
        error('measuring 0')
    end
    drSMA=2*norm(dxSMA,dySMA);%distance traveled to get to SMA surface and back
        xbounce=(2*dxSMA)+Pbounce(1);%new xbounce
        
        %determine the refractive indices && attenuation coefficients
        nt=SystemParam.n1;
        [ni,k]=medium_check(SystemParam,xbounce);%quick function to figure what medium we're in + the relevant proerties
        %     include a check %if the ray travels through the same medium the whole time
        %     if SystemParam.waterInterface==0 || ((Pbounce(1)<xbounce)&&(xbounce < SystemParam.waterstart)) || ((Pbounce(1)>xbounce)&&(xbounce > SystemParam.waterstart))
        kloss=exp(-k*(drSMA*10^-4));
        Vreturn=abs(Vbounce).*[sign(Vbounce(1)),-1*sign(Vbounce(2))];%return in new y direction
        %     else %if the ray switches into a new medium. need to calc losses,
        %     snells law, fresnell, new dir
        %         dx_int1=SystemParam.waterstart-Pbounce(1);
        %         dx_int2=xbounce-SystemParam.waterstart;
        %         if dx_int1<0%if we're traveling in the -x dir
        %             klosswater=
        %         if (Pbounce(1)+dxSMA) <SystemParam.waterstart%if the transition is in between the sma back to fiber travel path
        %         else
        %             Vi=Vbounce;%same initial direction
        %             dx_int=SystemParam.waterstart-(Pbounce(1)+dxSMA);
        %             dx_air=abs(dxSMA)+abs(dx_int);
        %             t_air=(abs(dxSMA)+abs(dx_int))/abs(Vbounce(1));
        %             dy_air=abs(t_air*Vbounce(2));
        %             dx_water=xbounce-SystemParam.waterstart;
        %             t_water=dx_water/abs(Vbounce(1));
        %             dy_water=abs(t_water*Vbounce(2));
        %         else
        
        I_bouncei=kloss*(1-metal_abs)*I_bounce;%input intensity into fresnel eq()
        I_bounceloss=(I_bounce-I_bouncei)+I_bounceloss;%loss amount
%%        [Tally,~,~] =  DifTrack(I_bounce,[I_bouncei,I_bounceloss],SystemParam,'BounceLoss w/in meas_SMA_check',1,Global_Index,Tally);

        %define the surface
        nhat=[0,-1];
        direction=sign(Vbounce(1));
        horz_surf=1;
        %snells law for angles and velocity
        %fresnell equation
        [theta_i,theta_t,theta_c,~,~,V_reflect,~] = Snells(Vreturn,nhat,ni,nt,horz_surf,direction);
        [I_reflect_sma,I_transmit_fib]=FresnelEq(I_bouncei,SystemParam,theta_i,theta_t,theta_c,ni,nt);
        %check if there's a lot of differences
%%        [Tally,~,~] =  DifTrack(I_bouncei,[I_reflect_sma,I_transmit_fib],SystemParam,'FresnellEq bounce w/in meas_SMA_check',1,Global_Index,Tally);
        
        %record the extra intensity coming back in
        I_extra_return(bounce)=I_transmit_fib;
%%        [Tally,~,~]=DifTrack(I_bounce,[I_extra_return(bounce),I_reflect_sma,I_bounceloss;],SystemParam,'bounceloop w/in meas_SMA_check',1,Global_Index,Tally);

        I_bounce=I_reflect_sma;
        Pbounce=[xbounce,P(2)];
        P_extra_return(bounce,:)=Pbounce;
        Vbounce=V_reflect;

        if I_bounce<SystemParam.I_min
            Tally.photon_min_count(iteration,h,aa)=Tally.photon_min_count(iteration,h,aa)+1;%counting all the times the model exits a function bc the direction is vertical
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
%%    [Tally,~,~]=DifTrack([I_transmit,I_reflect],[I_return,(Meas0-sum(meas.inten)),(House0-IT.housi),Cutoff0-IT.cutoffi;],SystemParam,'SMA fullbounce w/in meas_SMA_check',1,Global_Index,Tally);
%%    [Tally,~,~]=inOutTrack(I_in,(Meas0-sum(meas.inten)),SystemParam,'SMA fullbounce w/in meas_SMA_check',1,Global_Index,Tally);

    return
end
end
