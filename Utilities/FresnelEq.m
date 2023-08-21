function [I_reflect,I_transmit]=FresnelEq(I0,SystemParam,theta_i,theta_t,theta_c,ni,nt)
%using initial light intensity, incoming angle to the normal, transmitting
%angle from the normal, the critical angle, the refractive indices, then
%the distance traveled from last measurement to find the various
%reflection, transmission, abs. and backscat energies


    I1=I0;
    if I1==0
        I_reflect=0;
        I_transmit=0;
        return
    end
    
        % %transmitted light: fresnell's equations
    %     t_s = 2*ni*cos(theta_i)/((ni*cos(theta_i))+(nt*cos(theta_t))); %S polarized transmission
    %     t_p = 2*ni*cos(theta_i)/((ni*cos(theta_t))+(nt*cos(theta_i))); %P polarized transmission
    %     T_s=(nt*cos(theta_t)/(ni*cos(theta_i)))*t_s^2;%Transmittance =Transmitted/Incident power
    %     T_p=(nt*cos(theta_t)/(ni*cos(theta_i)))*t_p^2;

    %reflected light: fresnell's equations
        r_s=((ni*cos(theta_i))-(nt*cos(theta_t)))/((ni*cos(theta_i))+(nt*cos(theta_t)));
        r_p=((ni*cos(theta_t))-(nt*cos(theta_i)))/((ni*cos(theta_t))+(nt*cos(theta_i)));
        Rp=r_p^2; %
        Rs=r_s^2;


     %assuming equal amounts of s and p polarized light
    if isreal(theta_c)% if the cirtical angle is real (ie initial medium has greater RI than the transmitting medium
        %temporary using this. should remove or explain
        dX=SystemParam.cont_dx;
                    I_ray=I0*(1-10^-(SystemParam.alpha_r*dX/10));%light due to rayleigh loss
                    ray_ratio=I_ray/I0;
        theta_horiz=sign(theta_i)*pi/2-theta_i;%should change!
    RayleighCoeff=(1-10^(-SystemParam.alpha_r.*SystemParam.xlen/10));

    coeff=RayleighCoeff*2*(sin(theta_c)*cos(theta_c)*cos(2*theta_horiz)+3*theta_c);%(1+2*rand(1,1))

    forward_coeff=(1/2)*(sin(2*theta_c)*cos(2*theta_i)-6*theta_c+3*pi);%(1+2*rand(1,1));
    loss_coeff=(1/2)*((-cos(2*theta_horiz)*sin(2*theta_c))+6*theta_c+3*pi)/(3*pi);%(1+2*rand(1,1));
    %loss_coeff=(1/2)*((-cos(2*theta_i)*sin(2*theta_c))+6*theta_c+3*pi)/(3*pi);%(1+2*rand(1,1));

    forward_coeff=(3*pi)*(1-loss_coeff);
      back_coeff=forward_coeff/(3*pi);
      side_coeff=(loss_coeff-back_coeff)/(3*pi);
      scat_coeff=(ray_ratio)/(loss_coeff/(3*pi));
      coeff_norm=scat_coeff*side_coeff;
      I_backscat=(scat_coeff*back_coeff)*I0;
      
        if theta_i<theta_c || nt>=ni
        I_reflect = (((Rs+Rp)/2)*(I1*(1+coeff_norm)));%average of polarized light types reflected coeff* Intensity at interface
        I_transmit = ((I1*(1-coeff_norm))- I_reflect);%Transmitted/refracted power
        else
            I_transmit=coeff*I1;
            I_reflect=(1-coeff)*I1;

        end
    else %if the initial medium has lower RI than the transmitting medium then all of the transmitted light will go through
        I_reflect = (((Rs+Rp)/2)*(I1));
        I_transmit = (I1 - I_reflect);%Transmitted/refracted power
    end
end