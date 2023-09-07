function [I_reflect,I_transmit]=FresnelEq(I0,SystemParam,theta_i,theta_t,theta_c,theta_ih,ni1,nt1,horz_surf)
%using initial light intensity, incoming angle to the normal, transmitting
%angle from the normal, the critical angle, the refractive indices, then
%the distance traveled from last measurement to find the various
%reflection, transmission, abs. and backscat energies
%initial variables needed

%should be able to change if the mediums are lossy in the main function
lossy_f=SystemParam.lossy_f;
lossy_m=SystemParam.lossy_m;
%ang_wavenumber=(2*pi)/SystemParam.wavelength;%[rad/m]
wavenumber=1/(SystemParam.uv_wavelength*10^-9);%;[m-1]
c= 299792458;%[m/s] speed of light
ang_freq=wavenumber*2*pi*c;%[rad/s]
%attenuation coefficients in mediums. should be able to get this from a
%struct and/ or function
%change this
% alpha_uv=(0.9);%db/m friom https://www.content.molex.com/dxdam/literature/987650-8936.pdf
% alpha_r=(1.64*10^-3)*(850/275)^4;%db/m

%complex part of the refractive indexes
% a_uv=(alpha_uv/10)*log(10);%log(AUV);%(1/m)
% a_r=(alpha_r/10)*log(10);%;%-log(AR);%(1/m)
%summed attenuation in fiber
a_f=456.9589;%(a_uv+a_r);%*10^6;%[1/um*um/m]->[m-1]%for fiber
%attenuation coefficients in mediums
%relative permeability%TO CHANGE
mu_i=1;%;mu_i/mu_0;
mu_t=1;%mu_t/mu_0;

nf_i=a_f*c/(2*ang_freq);%imaginary portion of complex RI/ extinction coefficient
if any(imag(ni1))
    ni1=real(ni1);
    ni2=imag(ni1);
    lossy_i=lossy_m;
elseif ni1==real(SystemParam.n5)
    a_med=SystemParam.k*10^2;
    ni2=a_med*c/(2*ang_freq);
    lossy_i=lossy_m;
elseif ni1==real(SystemParam.n1)
    ni2=nf_i;
    lossy_i=lossy_f;
elseif ni1==real(SystemParam.n2)
    a_med=SystemParam.kair*10^2;
    ni2=a_med*c/(2*ang_freq);
    lossy_i=lossy_m;
else
    error('need more comprehensive coding of complex RI')
end
% % check=0;
if any(imag(nt1))
    nt1=real(nt1);
    nt2=imag(nt1);
    lossy_t=lossy_m;
% %     check=1;
elseif nt1==real(SystemParam.n1)
    nt2=nf_i;
    lossy_t=lossy_f;
% %     check=2;
elseif nt1==real(SystemParam.n5)
    a_med=SystemParam.k*10^2;
    nt2=a_med*c/(2*ang_freq);
    lossy_t=lossy_m;
% %     check=2;
elseif nt1==real(SystemParam.n2)
    a_med=SystemParam.kair*10^2;
    nt2=a_med*c/(2*ang_freq);
    lossy_t=lossy_m;
else
% %     check=3;

    error('need more comprehensive coding of complex RI')
end

%if the mediums are functionally the same, there's no change
if nt1==ni1
    I_reflect=0;
    I_transmit=I0;
    return
end
% % check2=0;
atten_ang_i=pi/2-theta_i;%incident attenution angle assumes homogenous plane wave
if lossy_i==0 && lossy_t==0
    %no imaginary part of the complex refraction index
    ni2=0;
    nt2=0;

    atten_ang_t=pi/2;%angle between attenuation coeff and surface is pi/2 nominally
% % check2=1;
elseif lossy_i==0 && lossy_t==1
    atten_ang_t=pi/2;%angle between attenuation coeff and surface is pi/2 nominally
       %no imaginary part of the complex refraction index
    ni2=0;
% %     check2=2;
elseif lossy_i==1 && lossy_t==0
    if any(nt2)
            atten_ang_t=(ni2/nt2).*sin(theta_i);%snells law attenuation angle
    else
        atten_ang_t=pi/2;
    end
     %no imaginary part of the complex refraction index
    nt2=0;
% %     check2=3;
elseif lossy_i==1 && lossy_t==1
% %     check2=4;
    if any(nt2)
        atten_ang_t=(ni2/nt2).*sin(theta_i);%snells law attenuation angle
    else
        atten_ang_t=pi/2;
    end
end
ei=(ni1^2-ni2^2)+(2*ni1*ni2)*1i;%complex relative permittivity of the incident
et=(nt1^2-nt2^2)+(2*nt1*nt2)*1i;%complex relative permittivity of the transmitted

    I1=I0;
    if I1==0
        I_reflect=0;
        I_transmit=0;
        return
    end
ni=ni1+1i*ni2;
nt=nt1+1i.*nt2;
ki=ni;%(ang_freq.*ni./c);
kt=nt;%(ang_freq.*nt./c);

kveci=[real(ki).*sin(theta_i)+(1i.*imag(ki).*cos(atten_ang_i));((real(ki).*cos(theta_i))+(1i.*imag(ki).*sin(atten_ang_i)))];
% kvecr=[kveci(1);((-1*real(ki).*cos(theta_i))+(1i.*imag(ki).*sin(atten_ang_i)))];
kvect=[real(ki).*sin(theta_i)+(1i.*imag(kt).*cos(atten_ang_t));((real(kt).*cos(theta_t))+(1i.*imag(kt).*sin(atten_ang_t)))];
    kiy=kveci(2);
%     kix=kveci(1);
% %     kry=kvecr(2);
%     krx=kvecr(1);
ev_fact=0;%evanescent wave factor ~= 1

if isreal(theta_c) && theta_i>=theta_c
    %kvect(1,i)=real(ki).*sin(theta_i)+(imag(kt).*cos(atten_ang_t));
    kvect(2)=1i.*(sqrt(((real(ki)*sin(theta_i))^2)-(real(kt).^2))+(imag(kt).*sin(atten_ang_t)));
    ev_fact=1;
end
kty=kvect(2);
ktx=kvect(1);
da_s=real((kty.*mu_i)./(kiy.*mu_t));
da_p=real((kty.*ei)./(kiy.*et));%(conj(kty).*ei.*conj(et))./(kiy.*et.*conj(et)));

rs=((mu_t.*kiy)-(mu_i.*kty))./((mu_t.*kiy)+(mu_i.*kty));
ts=(2.*mu_t.*kiy)./((mu_t.*kiy)+(mu_i.*kty));
rp=((ei.*kty)-(et.*kiy))./((et.*kiy)+(ei.*kty));
tp=(2.*et.*kiy)./((et.*kiy)+(ei.*kty));
        R_p=rp.^2; %
        R_s=rs.^2;
         T_s=da_s.*ts.^2;%Transmittance =Transmitted/Incident power
         T_p=da_p.*tp.^2;
        %check to make sure the real values work
        if real(R_s)>1
            R_s=1+imag(R_s);%should make a counter to track how often this happens....
        end
        if real(R_p)>1
            R_p=1+imag(R_p);
        end
if (real(R_s)<=0) && (real(R_p)<=0)
    R=0;
elseif (real(R_s)<=0) && (real(R_p)>0)
    R=real(R_p);
        if R>1
        R=1;
        %cshould count how many times this happens....
        end
elseif (real(R_p)<=0) && (real(R_s)>0)
    R=real(R_s);
        if R>1
        R=1;
        %cshould count how many times this happens....
        end
else
    R=(real(R_p)+real(R_s))/2;
    if R>1
        R=1;
        %cshould count how many times this happens....
    end
end

        I_reflect = (R*(I1));
        I_transmit = (I1 - I_reflect);%Transmitted/refracted power
        
        if I_reflect>I1 || I_reflect<0
            disp(I1)
            disp(I_reflect)
            disp(I_transmit)
            disp(R)
            disp(R_s)
            disp(R_p)
            disp(T_s)
            disp(T_p)
            error('I_reflect is oob')

        elseif I_transmit>I1 || I_transmit<0
            disp(I1)
            disp(I_reflect)
            disp(I_transmit)
            disp(R)
            disp(R_s)
            disp(R_p)
            disp(T_s)
            disp(T_p)
            error('I_transmit is oob')
       end
        e_decay=exp(-1*imag(kty));%exponential decay in y direciton. missing the 'y' exponent
        e_travel=exp(-1*1i.*ktx);%exponential decay in y direciton. missing the 'y' exponent

        E_wave_p=ev_fact.*tp.*I_reflect;
        E_wave_s=ev_fact.*ts.*I_reflect.*[0,0,1];%going in z direction
end
