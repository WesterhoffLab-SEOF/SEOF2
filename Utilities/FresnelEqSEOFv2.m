function [I_reflect,I_transmit]=FresnelEqSEOFv2(I0,theta_i,theta_t,theta_c,ni1,nt1)
%using initial light intensity, incoming angle to the normal, transmitting
%angle from the normal, the critical angle, the refractive indices, then
%the distance traveled from last measurement to find the various
%reflection, transmission, abs. and backscat energies
%initial variables needed

%Based on the lossy fresnell equations described by Shen & Kong Weber 1995, Weber 2014, and
%Oughtstun &Palombini, 2018

%syms n2p n2pp n1p n1pp theta_i psi_i psi_t mu_i mu_t ei et
n2p=real(nt1);
n2pp=imag(nt1);
n1p=real(ni1);
n1pp=imag(ni1);
%relative permeability% TO DO: CHANGE LATER FOR FULL E-M WAVE W CURRENT
mu_i=1;%;mu_i/mu_0;
mu_t=1;%mu_t/mu_0;

ni=n1p+1i*n1pp;
nt=n2p+1i.*n2pp;
ki=ni;%(ang_freq.*ni./c);
kt=nt;%(ang_freq.*nt./c);

%for a homogenous plane wave (according to oughsten & palombini 2018
psi_i=(pi/2)-theta_i;
if n1pp==0%if medium 1 is lossless, then psi_t=pi/2
    psi_t=pi/2;
else
    psi_t=acos((n1pp/n2pp).*cos(psi_i));%;%snells law attenuation angle using cosine bc sin(x)=cos(90-x)
end

%dielectric constant calculation
ei=(ni^2)/mu_i;%complex relative permittivity of the incident
et=(nt^2)/mu_t;%complex relative permittivity of the transmitted
I1=I0;
if I1==0%shouldn't happen, but exit in case
    I_reflect=0;
    I_transmit=0;
    return
end

if n1pp~=0 || n2pp~=0%if both mediums are lossy
    kveci=[n1p.*sin(theta_i)+(1i.*n1pp.*cos(psi_i));((n1p.*cos(theta_i))+(1i.*n1pp.*sin(psi_i)))];
    kvect=[n1p.*sin(theta_i)+(1i.*n2pp.*cos(psi_t));((n2p.*cos(theta_t))+(1i.*n2pp.*sin(psi_t)))];
    kiy=kveci(2);

    if n1pp~=0 && isreal(theta_c) && theta_i>=theta_c %in the supercritical condition, need to approximate stuff
        kvect(2)=1i.*(sqrt(((n1p*sin(theta_i))^2)-(n2p.^2))+(n2pp.*sin(psi_t)));
    end
    kty=kvect(2);

    rs=((mu_t.*kiy)-(mu_i.*kty))./((mu_t.*kiy)+(mu_i.*kty));
    rp=((et.*kiy)-(ei.*kty))./((et.*kiy)+(ei.*kty));

    Int_s=0;
    Int_p=0;
    R_p=norm(rp).^2;%abs(rp)^2;%.^2 %
    R_s=norm(rs).^2;%abs(rs)^2;%.^2
else %non complex , simpler fresnel code
    rp=(n1p*cos(theta_t)-n2p*cos(theta_i))./(n1p*cos(theta_t)+(n2p*cos(theta_i)));
    rs=(n1p*cos(theta_i)-n2p*cos(theta_t))./(n1p*cos(theta_i)+n2p*cos(theta_t));
    if theta_i>=theta_c && imag(theta_c)==0 %theta_c must be real and theta_i must exceed it
        R_p=1;
        R_s=1;
    else
        R_p=norm(rp)^2;%abs(rp)^2;%.^2 %
        R_s=norm(rs)^2;%abs(rs)^2;%.^2

    end

end

%average the reflection coefficient
R=(R_p+R_s)/2;
if R>1
    R=1;
end
I_reflect=R*I0;
I_transmit=I0-I_reflect;

