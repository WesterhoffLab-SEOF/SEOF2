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


iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);xx=Global_Index(4);yy=Global_Index(5);
%univeral parameters
n1 = SystemParam.n1;%fiber RI
% n2 = SystemParam.n2;%separation gap/air RI
% n5 = SystemParam.n5;%surrounding media  (water or air)RI
xlen = SystemParam.xLen;
intensityMin=SystemParam.intensityMin;
dtrav=SystemParam.contDx;%continuous transmission interval [uW]
%dxordr=SystemParam.dxordr;%select 0 if the logic is using dX, select 1 if the logic is using the whole travel distance

%define general surface
%TO DO: surface roughness update will require us to find nhat within the
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
while condition==1 && I>=intensityMin
    %record initial values for checking performance/accuracy
    MEAS0_win=sum(meas.inten);ABS0=IT.absorbi;BACK0=IT.backi;CUTOFF0=IT.cutoffi;HOUSE0=IT.housi;APPROX0=IT.approxi;TRANS0=IT.transi;B2H0=IT.backi;

    %%%%%%%make sure the ray isn't vertical
    Tally.whiletravcount(iteration,h,aa)=Tally.whiletravcount(iteration,h,aa)+1;
    if isnan(I)
        error('I is NaN in travleing')

    end
    xtoy=V(1)/V(2);
    if isreal(xtoy)&& xtoy==0%x vector is 0vector
        Meas_0=sum(meas.inten);Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;
        Tally.photon_vert_count(iteration,h,aa)=Tally.photon_vert_count(iteration,h,aa)+1;%counter goes up
        %%set up output variables
        Iend=0;
        Vend=[0,1*sign(V(2))];
        Pend=[P(1),bound(2)*sign(V(2))];

        %CALCULATE LOSSES
        %if it's within the part of the SMA connector, asummed that will
        %absorb everything
        if P(1)<=SystemParam.smaFlushLength && P(1)>=0
            %record the values
            %assume we lose all of the light into the housing, nothing is
            %measured
            IT.housi=I+IT.housi;
        else
            %record the measureable values. assume all light is side
            %emitted
            meas.points(meas.counter,:)=Pend;
            meas.inten(meas.counter,:)=I;
            meas.sum=meas.sum+I;
            meas.counter=meas.counter+1;
        end
        return%exit the function, go back to the main function
    end
    %parametric distance travel calc
    %check what the tbound in the y direction, and the tbound in the x
    %direction would be. if the tbound in the y direction is the shortest,
    %we use that then calculate a
    dx=dtrav*direction;t_std=dx/V(1);%standard travel t
    dxbound=bound(1)-(P(1));tb_x=dxbound/V(1);
    dybound=sign(V(2))*bound(2)-P(2);tb_y=dybound/V(2);
    tvec_i=[t_std,tb_x,tb_y];
    tvec=tvec_i;
    %check the parametric values are real and reasonable
    inf_ind=isinf(tvec);
    tvec(inf_ind)=10^10;%make it too large to be used
    %check the parametric values aren't negative (no such thing as negative
    %time) or 0
    tvec_neg_ind=tvec<=0;
    tvec(tvec_neg_ind)=10^10;%make it too large to be used

    if any(tvec_neg_ind) || isempty(tvec(tvec~=10^10))
        error('no directional travel')
    end
    [t_sort,t_idvec]=sort(tvec,'ascend');% choose the shortest path
    t=t_sort(1);%use the smallest t
    t_id=t_idvec(1);
    dx=t*V(1);dy=t*V(2);dr=norm([dx,dy]);
    x=dx+P(1);
    dx_bound_prev=dxbound;
    dy_bound_prev=dybound;
    tvec_prev=t_sort;
    %calculate the losses
    [nt,~,Tco]=medium_check(SystemParam,x);%quick function to figure what medium we're in + the relevant proerties
    ni=n1;
    switch t_id
        case 1 %travel forward without hitting any boundary
            Tally.case1count(iteration,h,aa)=Tally.case1count(iteration,h,aa)+1;

            P_in=[(P(1)+dx),P(2)+dy];%

            %loss calculation using dr
            %first need to do snell
            [theta_i,~,~,~,~,~,V_reflect,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit_calc,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,P_in,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            I_sideemit=Tco*I_sideemit_calc;
            %put side emitted light into the measure function
            if any(I_sideemit)%if there is side emitted light
                [I_return_m,~,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,V_reflect,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
                if isnan(I_return_m)
                    error('isNaN')
                end
            elseif any(I_sideemit_calc)%if there's not side emitted light but there could have been
                %then that light has been lost to the housing
                IT.housi=I_sideemit_calc+IT.housi;
            elseif isnan(I_sideemit)
                error('isNaN')
                %return
            end

            %create return values
            P_return=[P(1)+dx,P(2)+dy];
            V_return=V;
            V_return_sma=[0,0];
            I_return=I_remaining;
            if isnan(I_return)
                error('isNaN')
            end
        case 2 %travel forward until we hit the x boundary (tip of fiber)
            Tally.case2count(iteration,h,aa)=Tally.case2count(iteration,h,aa)+1;

            Pend=[(P(1)+dx),P(2)+dy];%
            P_in=Pend;
            %loss calculation using dr
            %first need to do snell
            [theta_i,~,~,~,~,~,V_reflect,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit_calc,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,P_in,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            I_sideemit=Tco*I_sideemit_calc ;
            %put side emitted light into the measure function
            if any(I_sideemit)%if there is side emitted light
                [I_return_m,~,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,V_reflect,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
                if isnan(I_return_m)
                    error('isNaN')
                end
            elseif any(I_sideemit_calc)%if there's not side emitted light but there could have been
                %then that light has been lost to the housing
                IT.housi=I_sideemit_calc+IT.housi;
            elseif ~any(I_sideemit_in)
                error('trying to measure 0 intensity in case 2')
            end
            %set up variables to exit the function
            Iend=I_remaining;
            if isnan(Iend)
                error('isNaN')
            end
            Vend=V;
            return
        case 3 %travel forward until we hit the y boundayr. need to calculate the changes
            Tally.case3count(iteration,h,aa)=Tally.case3count(iteration,h,aa)+1;
            P_in=[(P(1)+dx),sign(V(2))*bound(2)];%y boundary location

            %start with loss calculation using dr
            %start with snell. assign v_reflect as v_return
            [theta_i,theta_t,theta_c,~,~,V_transmit,V_reflect,direction] = Snells(V,nhat,ni,nt,horz_surf,direction);
            [I_remaining,I_sideemit_calc,I_sideemit_in,~,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V,I,P_in,dr,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index);
            I_sideemit=I_sideemit_calc*Tco;
            I_in=I_remaining;%set up remaining as input for the fresnell equation at the boundary
            %put cont loss side emitted light into the measure function. assum
            %nothing is returnng
            if any(I_sideemit)%if there is side emitted light
                [I_return_m,~,~,IT,Tally,meas]=meas_SMA_check(I_sideemit_in,I_sideemit,0,V_sideemit,V_reflect,P_sidemit,SystemParam,IT,Tally,meas,Global_Index);
                if isnan(I_return_m)
                    error('isNaN')
                end
            elseif any(I_sideemit_calc)%if there's not side emitted light but there could have been
                %then that light has been lost to the housing
                IT.housi=I_sideemit_calc+IT.housi;
            elseif ~any(I_sideemit_in)&& I_sideemit~=0
                error('trying to measure 0 intensity in case 3')
            end

            %now can calculate reflected and transmitted light using the
            %fresnell eq
            if I_in <0%=SystemParam.intensityMin
                error('I_in for fresnell case 3 <intensityMin')
            end
            [I_reflect_calc,I_transmit_calc]=FresnelEqSEOFv2(I_in,theta_i,theta_t,theta_c,ni,nt);

            I_transmit=I_transmit_calc*Tco;
            I_reflect=I_reflect_calc;
            I_ref_loss=0;
            %measure, if anything to transmit
            IT.housi=I_ref_loss+IT.housi;
            if  any(I_transmit)
                [I_return,P_return,V_return_sma,IT,Tally,meas]=meas_SMA_check(I_in,I_transmit,I_reflect,V_transmit,V_reflect,P_in,SystemParam,IT,Tally,meas,Global_Index);
                V_return=V_return_sma;

            elseif any(I_transmit_calc)%if there's not side emitted light but there could have been
                %then that light has been lost to the housing
                IT.housi=I_transmit_calc+IT.housi;
                I_return=I_reflect;
                P_return=P_in;
                V_return=V_reflect;
                V_return_sma=[0,0];

            elseif any(I_in)
                I_return=I_reflect;
                P_return=P_in;
                V_return=V_reflect;
                V_return_sma=[0,0];
            else
                error('trying to measure 0 intensity in case 3')
            end
        otherwise
            error('no good case')
    end
    %record starting values for next while loop
    P=P_return;
    x=P_return(1);
    V=V_return;
    I=I_return;

    if isnan(I)
        error('isNaN')
    end

    if I<SystemParam.intensityMin
        Tally.minPhotons_count(iteration,h,aa)=Tally.minPhotons_count(iteration,h,aa)+1;%counting all the times the model exits a function bc the direction is vertical
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

%if while loop terminates and end values aren't assigned
Iend=I;
if isnan(Iend)
    error('unexpected direction')
end
Pend=P;
Vend=V;

end

