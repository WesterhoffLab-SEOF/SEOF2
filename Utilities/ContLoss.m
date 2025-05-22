function [I_remaining,I_sideemit,I_sideemit_in,I_losstot,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V_in,I_in,P,dtrav,bound,ni1,nt1,theta_i,SystemParam,Tally,IT,Global_Index)
%function for applying the continuous losses, including side scattering and absorption losses, at every dx
%[I_remaining,I_sideemit,P_sidemit,IT,Tally]=ContLoss(V_in,I_in,P,[dx,dy],bound,ni,nt,nhat,theta_i,theta_c,SystemParam,Tally,IT,Global_Index]
%starting metrics for tracking possible differences
Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;
%CURRENTLY ONLY APPLICABLE FOR HORIZONTAL SURFACES
horz_surf=1;
%make sure this is real
ni=real(ni1);
nt=real(nt1);
ni_c=ni1;
nt_c=nt1;
dt=norm(dtrav);
alpha=SystemParam.alpha;%db/um
r_coeff=SystemParam.rayScatterCoeff;%coefficient of losses attributed to rayleigh
losses=(I_in*(1-10^(-alpha*dt/10)));
Absorb_in=losses*(1-r_coeff);

%calculate remaining light after absorption. wait to remove more 
I_rem1=I_in-Absorb_in;
%record absorbed light.
IT.absorbi=IT.absorbi+Absorb_in;

%start considering rayleigh scattered light
%need to find appropriate indices to get the pre calculated scattering
%coefficients
deg=abs(rad2deg(theta_i));
deg_index=round(SystemParam.angleDiv*deg);%index for the angle
%find the RI index 
switch ni
    case real(SystemParam.n1)
        ni_index=1;
    case real(SystemParam.n2)
        ni_index=2;
    case real(SystemParam.n3)
        ni_index=3;
    case real(SystemParam.n4)
        ni_index=4;
    case real(SystemParam.nWater)
        ni_index=5;
    case real(SystemParam.n_metal)
        ni_index=6;
end

switch nt
    case real(SystemParam.n1)
        nt_index=1;
    case real(SystemParam.n2)
        nt_index=2;
    case real(SystemParam.n3)
        nt_index=3;
    case real(SystemParam.n4)
        nt_index=4;
    case real(SystemParam.nWater)
        nt_index=5;
    case real(SystemParam.nMetal)
        nt_index=6;
end
%import the pre calculated scattering coefficients
scat_coeffs=cell2mat(SystemParam.scatterCoeff(ni_index,nt_index,deg_index));
%assign the coefficients to appropriate var
Forward=scat_coeffs(1);Backward=scat_coeffs(2);Sideup=scat_coeffs(3);Sidedown=scat_coeffs(4);Rayloss=scat_coeffs(5);
%calculate amount of light scattered in each direction
%assumed lost Rayleigh light
RCloss_in=losses*r_coeff;%amount 'lost (ie not going forward) along the travel distance
I_scat_tot=RCloss_in.*Rayloss;%take the ratio for hwo much total light is actually scattered
I_forward=Forward*I_scat_tot;%add this back to the returning
I_backscat=Backward*I_scat_tot;%this is the main light reported to be lost to general scatter
I_scatup=I_scat_tot*Sideup;
I_scatdown=I_scat_tot*Sidedown;
I_rem2=I_rem1-I_scat_tot;%remaining light after the scattered amount is removed
%determine the approx distance the scatter light travels in
% disp(Forward+Backward+Sideup+Sidedown)
ydir=sign(V_in(2));
if ydir==0
    ydir=1;
end
dy_scatup=ydir*bound(2)-P(2);
dy_scatdown=-1*ydir*bound(2)-P(2);
Y_down=-1*ydir*bound(2);
Y_up=ydir*bound(2);
%assume all scattered light is in the vertical direction
theta_i=0;theta_t=0;theta_c=asin(nt/ni);theta_ih=pi/2;
%calculate the amount of scattered light lost on the way to the surface
%boundary. include losses for both absorb and rayleigh
I_scat_absorbedup=(I_scatup*(1-10^(-alpha*dy_scatup/10)));
I_scat_absorbeddown=(I_scatdown*(1-10^(-alpha*dy_scatdown/10)));
IT.absorbi=I_scat_absorbedup+I_scat_absorbeddown+IT.absorbi;%record as absorbed light
%remove the loss from up and down
I_scat_inup=I_scatup*(10^(-alpha*dy_scatup/10));
I_scat_indown=I_scatdown*(10^(-alpha*dy_scatdown/10));
%do fresnell eq to see what transmits
if any(I_scat_inup)% if the value isn't 0
    [I_scatup_reflect,I_scatup_transmit]=FresnelEq(I_scat_inup,SystemParam,theta_i,theta_t,theta_c,theta_ih,ni_c,nt_c,horz_surf);
else
    I_scatup_reflect=0;I_scatup_transmit=0;
end
if any(I_scat_indown)
    [I_scatdown_reflect,I_scatdown_transmit]=FresnelEq(I_scat_indown,SystemParam,theta_i,theta_t,theta_ih,theta_c,ni_c,nt_c,horz_surf);
else
    I_scatdown_reflect=0;I_scatdown_transmit=0;
end
%track differences here if there's an issue
% % [Tally,~,~] =  DifTrack(I_scatup,[I_scat_absorbedup,I_scatup_reflect,I_scatup_transmit],SystemParam,'Intensity for Scatup in ContLoss',1,Global_Index,Tally);
% % [Tally,~,~] =  DifTrack(I_scatdown,[I_scat_absorbeddown,I_scatdown_reflect,I_scatdown_transmit],SystemParam,'Intensity for Scatdown in ContLoss',1,Global_Index,Tally);


%won't continue to follow the reflected scatter light, record as an
%approximation
IT.approxi=I_scatup_reflect+I_scatdown_reflect+IT.approxi;
%record the backscattered light as the backscatter loss
IT.backi=I_backscat+IT.backi;

%tally the losses to get I remaining
I_scatterloss=I_backscat+I_scatup_reflect+I_scatdown_reflect+I_scat_absorbedup+I_scat_absorbeddown;%things backscattered or assumed to be very little and not worth following
I_scatteruse=I_scatup_transmit+I_scatdown_transmit;
I_scatteradd=I_forward;

%create variables for exiting the fiber
I_remaining=I_rem2+I_scatteradd;
%
I_sideemit0=[I_scatup_transmit;I_scatdown_transmit];
P_sidemit0=[P(1),Y_up;P(1),Y_down];
V_sideemit0=[0,ydir;0,-1*ydir];
[~,ind]=max(I_sideemit0);
%record up + down together to save time. pick higher value side emission to
%use for P and V data
I_sideemit=sum(I_sideemit0);%
I_sideemit_in=I_scat_inup+I_scat_indown;%
P_sidemit=P_sidemit0(ind,:);
V_sideemit=V_sideemit0(ind,:);
I_losstot=I_scatterloss+Absorb_in;
if isnan(I_remaining)
    error('nan in conttrav')
elseif ~any(I_remaining)
    error('I remaining is empty')
end
% %  Pow_Out=[I_remaining,I_sideemit,(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];
% % [Tally,~,~] =  DifTrack(I_in,Pow_Out,SystemParam,'ContLoss total',1,Global_Index,Tally);
% % [Tally,~,~] =  inOutTrack(I_in,I_sideemit,SystemParam,'ContLoss total',1,Global_Index,Tally);
if I_remaining<0
    error('Negative I remaining in conttrav')
end
