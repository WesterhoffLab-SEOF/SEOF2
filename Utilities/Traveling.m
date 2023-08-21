function [Iend,Vend,Pend,direction,meas,IT,Tally] = Traveling(I,V,P,direction,SystemParam,Bounds,Global_Index,Tally,meas)
%ray traveling along the fiber length
%%%%%initializing variables before the while loop
%temporary storage zero values
IT=travel_storage;%temporary storage of light losses within the travel function
IT.absorbi=0;
IT.backi=0;
IT.cutoffi=0;
IT.housi=0;
IT.transi=0;
IT.b2hi=0;
IT.approxi=0;

MEAS0=sum(meas.inten);
Istart=I;
SystemParam.SMA=1;%is there a SMA connector used to couple thhe fiber to the LED? 1 for yes, 0 for no
SystemParam.SMA_flushlength=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long
SystemParam.SMA_totallength=2.5*10^4;%the total length of the SMA connector is 2.5 cm long
SystemParam.SMA_diam=2.5*10^3;%diameter of SMA connector [mm to um
iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);
%univeral parameters
n1 = SystemParam.n1;%fiber RI
% n2 = SystemParam.n2;%separation gap/air RI
% n5 = SystemParam.n5;%surrounding media  (water or air)RI
xlen = SystemParam.xlen;
I_min=SystemParam.I_min;
dtrav=SystemParam.cont_dx;%continuous transmission interval [uW]
%dxordr=SystemParam.dxordr;%select 0 if the logic is using dX, select 1 if the logic is using the whole travel distance

%define general surface
%nhat should be moved to be found within the while loop when surface
%roughness included
    nhat=[0,-1];%surface normal
    horz_surf=1;%predominantly considering a horizontal surface
x=P(1);
%switch case statements save time bc matlab will just jump directly to the
%applicable statement. can't use it for complicated statements
%check conditions
switch direction
    case -1
        condition=(x>0);
        bound=Bounds.Pf0;%begining of fiber boundary
    case 1
        condition=(x<xlen);
        bound=Bounds.Pfe;%end of fiber boundary
    otherwise
        error('unexpected direction')
end
while condition && I>=I_min
    %record initial values for checking performance/accuracy
    MEAS0_win=sum(meas.inten);ABS0=IT.absorbi;BACK0=IT.backi;CUTOFF0=IT.cutoffi;HOUSE0=IT.housi;APPROX0=IT.approxi;TRANS0=IT.transi;B2H0=IT.backi;

    %%%%%%%make sure the ray isn't vertical
    Tally.whiletravcount(iteration,h,aa)=Tally.whiletravcount(iteration,h,aa)+1;
    if isnan(I)
        %disp(Tally.whiletravcount(iteration,h,aa))
        %disp(IT.absorbi)
    %disp(IT.backi)
    %disp(IT.cutoffi)
    %disp(IT.housi)
    %disp(V)
    %disp(P)
    error('I is bad in travleing')
    
    end
    xtoy=V(1)/V(2);
    if isreal(xtoy)&& xtoy==0%x vector is 0.1% or less of the y vector
        Meas_0=sum(meas.inten);Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;

        Tally.photon_vert_count(iteration,h,aa)=Tally.photon_vert_count(iteration,h,aa)+1;%counter goes up
        %%set up output variables
        Iend=0;
        %theta_end=pi/2;
        Vend=[0,1*sign(V(2))];
        Pend=[P(1),bound(2)*sign(V(2))];

        %CALCULATE LOSSES
        %if it's within the part of the SMA connector that will absorb all
        %of it. (But hey. maybe it'll scatter some stuff? unclear)
        if P(1)<=SystemParam.SMA_flushlength && P(1)>=0
            %record the values
            %assume we lose all of the light into the housing, nothing is
            %measured
            IT.housi=I+IT.housi;
        else
            %record the measureable values. assume all light is side
            %emitted
            meas.points(meas.counter,:)=Pend;
            meas.inten(meas.counter,:)=I;
            meas.counter=meas.counter+1;
        end
        Pow_Out=[Iend,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
        [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'vertical cond  w/in trav',1,Global_Index,Tally);
        [Tally,~,~] =  inOutTrack(I,(sum(meas.inten)-Meas_0),SystemParam,'vertical cond  w/in trav',1,Global_Index,Tally);

        return%exit the function, go back to the main function
    end
    %parametric distance travel calc
    %check what the tbound in the y direction, and the tbound in the x
    %direction would be. if the tbound in the y direction is the shortest,
    %we use that then calculate a 
    dx=dtrav*direction;t_std=dx/V(1);%standard travel t
    dxbound=bound(1)-(P(1));tb_x=dxbound/V(1);
    dybound=sign(V(2))*bound(2)-P(2);tb_y=dybound/V(2);
    tvec=[t_std,tb_x,tb_y];
    %check the parametric values are real and reasonable
    inf_ind=isinf(tvec);
    tvec(inf_ind)=10^10;%make it too large to be used
    if isempty(tvec)
        error('no directional travel')
    end
    [t_sort,t_idvec]=sort(tvec,'ascend');% choose the shortest path
    t=t_sort(1);%use the smallest t
    t_id=t_idvec(1);
    dx=t*V(1);dy=t*V(2);dr=norm([dx,dy]);
    %calculate the losses
    [nt,~]=medium_check(SystemParam,x);%quick function to figure what medium we're in + the relevant proerties
    ni=n1;
    switch t_id
        case 1 %travel forward without hitting any boundary
            Tally.case1count(iteration,h,aa)=Tally.case1count(iteration,h,aa)+1;
            Meas_0=sum(meas.inten);Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;

            P_in=[(P(1)+dx),P(2)+dy];%

            %loss calculation using dr
            %first need to do snell
            [theta_i,~,~,~,~,~,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,P_in,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            Pow_Out=[I_remaining,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'ContLoss w/in case 1',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I,I_sideemit,SystemParam,'ContLoss w/in case 1',1,Global_Index,Tally);

            %put side emitted light into the measure function
            if any(I_sideemit_in) && any(I_sideemit)
                meas_0=sum(meas.inten);cutoff0=IT.cutoffi;abs0=IT.absorbi;house0=IT.housi;approx0=IT.approxi;back0=IT.backi;trans0=IT.transi;b2h0=IT.b2hi;
            [I_return_meas,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
            Pow_Out=[I_return_meas,(sum(meas.inten)-meas_0),(IT.cutoffi-cutoff0),(IT.absorbi-abs0),(IT.housi-house0),(IT.approxi-approx0),(IT.backi-back0),(IT.transi-trans0),(IT.b2hi-b2h0)];
            [Tally,~,~] =  DifTrack(I_sideemit_in,Pow_Out,SystemParam,'meas_SMA_check  w/in case 1',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I_sideemit_in,sum(meas.inten)-Meas_0,SystemParam,'meas_SMA_check  w/in case 1',1,Global_Index,Tally);

            else
                error('trying to measure 0 intensity in case 1')
            end
    
            %create return values
            P_return=[P(1)+dx,P(2)+dy];
            V_return=V;
            I_return=I_remaining;
            Pow_Out=[I_return,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'case 1 total',1,Global_Index,Tally);
            [Tally,~,~] =  DifTrack(I,sum(meas.inten)-Meas_0,SystemParam,'case 1 total',1,Global_Index,Tally);

        case 2 %travel forward until we hit the x boundary
            Tally.case2count(iteration,h,aa)=Tally.case2count(iteration,h,aa)+1;
            Meas_0=sum(meas.inten);Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;

            Pend=[(P(1)+dx),P(2)+dy];%
            %loss calculation using dr
            %first need to do snell
            [theta_i,~,~,~,~,~,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,Pend,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            Pow_Out=[I_remaining,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'ContLoss w/in case 2',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I,I_sideemit,SystemParam,'ContLoss w/in case 2',1,Global_Index,Tally);

            %put side emitted light into the measure function
            if any(I_sideemit_in) && any(I_sideemit)
                meas_0=sum(meas.inten);cutoff0=IT.cutoffi;abs0=IT.absorbi;house0=IT.housi;approx0=IT.approxi;back0=IT.backi;trans0=IT.transi;b2h0=IT.b2hi;
                [I_return_meas,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
                Pow_Out=[I_return_meas,(sum(meas.inten)-meas_0),(IT.cutoffi-cutoff0),(IT.absorbi-abs0),(IT.housi-house0),(IT.approxi-approx0),(IT.backi-back0),(IT.transi-trans0),(IT.b2hi-b2h0)];
                [Tally,~,~] =  DifTrack(I_sideemit_in,Pow_Out,SystemParam,'meas_SMA_check  w/in case 2',1,Global_Index,Tally);
                [Tally,~,~] =  inOutTrack(I,sum(meas.inten)-Meas_0,SystemParam,'meas_SMA_check  w/in case 2',1,Global_Index,Tally);

            else                
                error('trying to measure 0 intensity')
            end
            %set up variables to exit the function
            Iend=I_remaining;
            Vend=V;
            Pow_Out=[Iend,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'case 2 total',1,Global_Index,Tally);
            [Tally,~,~] =  DifTrack(I,(sum(meas.inten)-Meas_0),SystemParam,'case 2 total',1,Global_Index,Tally);

            return
        case 3 %travel forward until we hit the y boundayr. need to calculate the changes
            Tally.case3count(iteration,h,aa)=Tally.case3count(iteration,h,aa)+1;
            Meas_0=sum(meas.inten);Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;


            P_in=[(P(1)+dx),sign(V(2))*bound(2)];%y boundary location

            %start with loss calculation using dr
            %start with snell. assign v_reflect as v_return
            [theta_i,theta_t,theta_c,~,V_transmit,V_reflect,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,P_in,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            Pow_Out=[I_remaining,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'ContLoss w/in case 3',1,Global_Index,Tally);
             [Tally,~,~] =  DifTrack(I,I_sideemit,SystemParam,'ContLoss w/in case 3',1,Global_Index,Tally);

             V_return=V_reflect;
             I_in=I_remaining;%set up remaining as input for the fresnell equation at the boundary
            %put cont loss side emitted light into the measure function. assum
            %nothing is returnng
            if any(I_sideemit_in) && any(I_sideemit)
                meas_0=sum(meas.inten);cutoff0=IT.cutoffi;abs0=IT.absorbi;house0=IT.housi;approx0=IT.approxi;back0=IT.backi;trans0=IT.transi;b2h0=IT.b2hi;

            [I_return_meas,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
            Pow_Out=[I_return_meas,(sum(meas.inten)-meas_0),(IT.cutoffi-cutoff0),(IT.absorbi-abs0),(IT.housi-house0),(IT.approxi-approx0),(IT.backi-back0),(IT.transi-trans0),(IT.b2hi-b2h0)];
            [Tally,~,~] =  DifTrack(I_sideemit_in,Pow_Out,SystemParam,'meas_SMA_check  w/in case 3,1',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I,(sum(meas.inten)-meas_0),SystemParam,'meas_SMA_check  w/in case 3,1',1,Global_Index,Tally);

            else
                %disp(I)
                %disp(I_sideemit)
                %disp(I_sideemit_in)
                %disp(I_losstot)
                %disp(P_sidemit)
                %disp(V_sideemit)
                
                error('trying to measure 0 intensity')
            end
            %now can calculate reflected and transmitted light using the
            %fresnell eq
            [I_reflect,I_transmit]=FresnelEq(I_in,SystemParam,theta_i,theta_t,theta_c,ni,nt);
            [Tally,~,~] =  DifTrack(I_in,[I_reflect,I_transmit],SystemParam,'FresnelEq w/in case 3',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I_in,I_transmit,SystemParam,'FresnelEq w/in case 3',1,Global_Index,Tally);

            %measure, if anything to transmit
            if any(I_in) || any(I_transmit) || any(I_reflect)
                %starting values
                meas_0=sum(meas.inten);cutoff0=IT.cutoffi;abs0=IT.absorbi;house0=IT.housi;approx0=IT.approxi;back0=IT.backi;trans0=IT.transi;b2h0=IT.b2hi;

            [I_return,P_return,IT,Tally,meas]=meas_SMA_check(I_in,I_transmit,I_reflect,V_transmit,P_in,SystemParam,IT,Tally,meas,Global_Index);
            %ending values
            Pow_Out=[I_return,(sum(meas.inten)-meas_0),(IT.cutoffi-cutoff0),(IT.absorbi-abs0),(IT.housi-house0),(IT.approxi-approx0),(IT.backi-back0),(IT.transi-trans0),(IT.b2hi-b2h0)];
                        %check if big dif
            [Tally,~,~] =  DifTrack(I_in,Pow_Out,SystemParam,'meas_SMA_check w/in case 3,2',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I,(sum(meas.inten)-Meas_0),SystemParam,'meas_SMA_check w/in case 3,2',1,Global_Index,Tally);

            else

                
                error('trying to measure 0 intensity')
            end
            Pow_Out=[I_return,(sum(meas.inten)-Meas_0),(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'case 3 total',1,Global_Index,Tally);
            
        otherwise
            error('no good case')
    end
  %check if big differences w/in while loop
            Pow_Out=[I_return,(sum(meas.inten)-MEAS0_win),(IT.cutoffi-CUTOFF0),(IT.absorbi-ABS0),(IT.housi-HOUSE0),(IT.approxi-APPROX0),(IT.backi-BACK0),(IT.transi-TRANS0),(IT.b2hi-B2H0)];

            [Tally,~,~] =  DifTrack(I,Pow_Out,SystemParam,'while loop w/in traveling',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(I,(sum(meas.inten)-MEAS0_win),SystemParam,'while loop w/in traveling',1,Global_Index,Tally);


  %record starting values for next while loop
  P=P_return;
  x=P_return(1);
  V=V_return;
  I=I_return;
 
    if I<SystemParam.I_min             
        Tally.photon_min_count(iteration,h,aa)=Tally.photon_min_count(iteration,h,aa)+1;%counting all the times the model exits a function bc the direction is vertical
        break
    end
   %check conditions
   switch direction
       case -1
           condition=(x>0);
           bound=Bounds.Pf0;%begining of fiber boundary
       case 1
           condition=(x<xlen);
           bound=Bounds.Pfe;%end of fiber boundary
       otherwise
           error('unexpected direction')
   end
    
end
  %check if big differences after while loop is over
            Pow_Out=[I_return,(sum(meas.inten)-MEAS0),IT.cutoffi,IT.housi,IT.approxi,IT.backi,IT.transi,IT.b2hi];
            [Tally,~,~] =  DifTrack(Istart,Pow_Out,SystemParam,'while loop total w/in traveling',1,Global_Index,Tally);
            [Tally,~,~] =  inOutTrack(Istart,(sum(meas.inten)-MEAS0),SystemParam,'while loop total w/in traveling',1,Global_Index,Tally);

%if while loop terminates and end values aren't assigned
Iend=I;
Pend=P;
Vend=V;

end

