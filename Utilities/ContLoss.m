function [I_remaining,I_sideemit,I_sideemit_in,I_losstot,V_sideemit,P_sidemit,IT,Tally]=ContLoss(V_in,I_in,P,dtrav,bound,ni,nt,theta_i,SystemParam,Tally,IT,Global_Index)
%[I_remaining,I_sideemit,P_sidemit,IT,Tally]=ContLoss(V_in,I_in,P,[dx,dy],bound,ni,nt,nhat,theta_i,theta_c,SystemParam,Tally,IT,Global_Index]
%starting metrics for tracking possible differences
Cutoff0=IT.cutoffi;Abs0=IT.absorbi;House0=IT.housi;Approx0=IT.approxi;Back0=IT.backi;Trans0=IT.transi;B2h0=IT.b2hi;


%find I remaining of the main ray afte scattered and absorption losses are
%applied

dt=norm(dtrav);
alpha_r=SystemParam.alpha_r;%db/um
alpha_uv=SystemParam.alpha_uv;%db/um
Absorb_in=(I_in*(1-10^(-alpha_uv*dt/10)));

%calculate remaining light after absorption. wait to remove more 
I_rem1=I_in-Absorb_in;
%record absorbed light.
IT.absorbi=IT.absorbi+Absorb_in;

%start considering rayleigh scattered light
%need to find appropriate indices to get the pre calculated scattering
%coefficients
deg=abs(rad2deg(theta_i));
deg_index=round(SystemParam.ang_div*deg);%index for the angle
%find the RI index 
switch ni
    case 1
        ni_index=2;
    case 1.333
        ni_index=5;
    case 1.353
        ni_index=4;
    case 1.5
        ni_index=1;
end

switch nt
    case 1
        nt_index=2;
    case 1.333
        nt_index=5;
    case 1.353
        nt_index=4;
    case 1.5
        nt_index=1;
end
%import the pre calculated scattering coefficients
scat_coeffs=cell2mat(SystemParam.scatter_coeff(ni_index,nt_index,deg_index));
%assign the coefficients to appropriate var
Forward=scat_coeffs(1);Backward=scat_coeffs(2);Sideup=scat_coeffs(3);Sidedown=scat_coeffs(4);Rayloss=scat_coeffs(5);
%calculate amount of light scattered in each direction
%assumed lost Rayleigh light
RCloss_in=I_in*(1-10^(-alpha_r*dt/10));%amount 'lost (ie not going forward) along the travel distance
I_scat_tot=RCloss_in.*Rayloss;%take the ratio for hwo much total light is actually scattered
I_forward=Forward*I_scat_tot;%add this back to the returning
I_backscat=Backward*I_scat_tot;%this is the main light reported to be lost to general scatter
I_scatup=I_scat_tot*Sideup;
I_scatdown=I_scat_tot*Sidedown;
I_rem2=I_rem1-I_scat_tot;%remaining light after the scattered amount is removed
%determine the approx distance the scatter light travels in
ydir=sign(V_in(2));
if ydir==0
    ydir=1;
end
dy_scatup=ydir*bound(2)-P(2);
dy_scatdown=-1*ydir*bound(2)-P(2);
Y_down=-1*ydir*bound(2);
Y_up=ydir*bound(2);
%assume all scattered light is in the vertical direction
theta_i=0;theta_t=0;theta_c=asin(nt/ni);
%calculate the amount of scattered light absorbed on the way to the surface
%boundary
I_scat_absorbedup=(I_scatup*(1-10^(-alpha_uv*dy_scatup/10)));
I_scat_absorbeddown=(I_scatdown*(1-10^(-alpha_uv*dy_scatdown/10)));
IT.absorbi=I_scat_absorbedup+I_scat_absorbeddown+IT.absorbi;%record as absorbed light
%remove the absorbed from up and down
I_scat_inup=I_scatup*(10^(-alpha_uv*dy_scatup/10));
I_scat_indown=(I_scatdown*(10^(-alpha_uv*dy_scatdown/10)));
%do fresnell eq to see what transmits
if any(I_scat_inup)% if the value isn't 0
    [I_scatup_reflect,I_scatup_transmit]=FresnelEq(I_scat_inup,SystemParam,theta_i,theta_t,theta_c,ni,nt);
else
    I_scatup_reflect=0;I_scatup_transmit=0;
end
if any(I_scat_indown)
    [I_scatdown_reflect,I_scatdown_transmit]=FresnelEq(I_scat_indown,SystemParam,theta_i,theta_t,theta_c,ni,nt);
else
    I_scatdown_reflect=0;I_scatdown_transmit=0;
end
%track differences here if there's an issue
[Tally,~,~] =  DifTrack(I_scatup,[I_scat_absorbedup,I_scatup_reflect,I_scatup_transmit],SystemParam,'Intensity for Scatup in ContLoss',1,Global_Index,Tally);
[Tally,~,~] =  DifTrack(I_scatdown,[I_scat_absorbeddown,I_scatdown_reflect,I_scatdown_transmit],SystemParam,'Intensity for Scatdown in ContLoss',1,Global_Index,Tally);


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
end
Pow_Out=[I_remaining,I_sideemit,(IT.cutoffi-Cutoff0),(IT.absorbi-Abs0),(IT.housi-House0),(IT.approxi-Approx0),(IT.backi-Back0),(IT.transi-Trans0),(IT.b2hi-B2h0)];

[Tally,~,~] =  DifTrack(I_in,Pow_Out,SystemParam,'ContLoss total',1,Global_Index,Tally);
[Tally,~,~] =  inOutTrack(I_in,I_sideemit,SystemParam,'ContLoss total',1,Global_Index,Tally);

end
