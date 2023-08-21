% function [I_scat,Theta_scat,V_scat,I_trans,direction,b2h]= EndReflect(I_in,V_in,SystemParam,directioni)
function [I_scat,Theta_scat,I_trans,direction,b2h,bounce_num,Tally,meas]= EndReflect(I_in,V_in,P_in,SystemParam,directioni,bounce_num,Global_Index,Tally,meas)
%takes incoming intensity,velocity,direction,number ofbounces. outputs the
%transmitted light intensity, scattered light intensity +direction+angle.
%Keeps track of the 
%define the surface interface
horz_surf=0;%vertical surface
vert_surf=1;%vertical surface
nhat=[-1,0];%assume it's perfectly flat (include surface roughness)
ni=SystemParam.n1;%incoming media is the fiber

if directioni==1%hitting the end of the fiber
    nt=SystemParam.n5;%transmitting media is whatever the fiber is submerged in

else%hitting the coupling face of the fiber
    nt=SystemParam.n2;%RI of "gap"/air %air index
end
%global index
 Meas0=sum(meas.inten);
iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);%indexes
%make sure the ray isn't vertical
xtoy=V_in(1)/V_in(2);
if isreal(xtoy)&& xtoy==0%if the ray is entirely vertical, it only side emits
    Tally.photon_vert_count(iteration,h,aa)=Tally.photon_vert_count(iteration,h,aa)+1;%counter goes up
    %record the values
    meas.points(meas.counter,:)=P_in;
    meas.inten(meas.counter,:)=I_in;
    meas.counter=meas.counter+1;
    meas.sum=I_in+meas.sum;
    %set up output variables
    I_trans=0;
    I_scat=0;
    Theta_scat=0;
    direction=directioni;
    b2h=0;
    [Tally,~,~] = DifTrack(I_in,[I_scat,I_trans,(sum(meas.inten)-Meas0)],SystemParam,'Vertical cond w/in EndReflect',1,Global_Index,Tally);
        [Tally,~,~] = inOutTrack(I_in,(sum(meas.inten)-Meas0),SystemParam,'Vertical cond w/in EndReflect',1,Global_Index,Tally);

    return%exit the function, go back to the main function
end
[theta_i,theta_t,theta_c,theta_r,~,~,direction] = Snells(V_in,nhat,ni,nt,horz_surf,directioni);
%direction here should flip from the original value due to the x directions
if directioni*direction~=-1
    error("direction didn't flip inside the end reflect value")
end
%of the incoming and leaving velocities being opposite
[I_reflect,Itrans]=FresnelEq(I_in,SystemParam,theta_i,theta_t,theta_c,ni,nt);
%check there isn't a large difference after using the fresnell eq. keep
%track in Tally struct
[Tally,~,~] =  DifTrack(I_in,[I_reflect,Itrans],SystemParam,'FresnelEq w/in EndReflect',1,Global_Index,Tally);
[Tally,~,~] =  inOutTrack(I_in,Itrans,SystemParam,'FresnelEq w/in EndReflect',1,Global_Index,Tally);

if directioni==1%hitting the end of the fiber
    I_trans=Itrans;%transmitted light will go into the medium
    b2h=0;
else%hitting the coupling face of the fiber
    b2h=Itrans;%transmitted light will fo back to the houing
    I_trans=0;
end
%find an appropriately sized scatter cone
[I_scatter,Theta_scat]=scatter_cone(SystemParam,theta_r,vert_surf);
I_scat=I_scatter.*I_reflect;
%check there isn't a large difference after scattering. keep
%track in Tally struct
[Tally,~,~] =  DifTrack(I_reflect,I_scat,SystemParam,'scatter_cone inside EndReflect',1,Global_Index,Tally);

%update number of bounces
bounce_num=bounce_num+1;

[Tally,~,~] =  DifTrack(I_in,[(sum(meas.inten)-Meas0),I_scat,I_trans,b2h],SystemParam,'EndReflect total',1,Global_Index,Tally);

end
