function [Ray_X,Ray_Y,alpha_ang,beta_ang,Intensity] = OutputLED3D(SystemParam)
%initial variable set up using System Param Struct
LED_height=SystemParam.ledDiameter*10^3;%(mm>um)
LED_width=SystemParam.ledDiameter*10^3;%(mm>um)
ledAngle=SystemParam.ledAngle;
ledHalfAngle=SystemParam.ledHalfAngle;
P_init=SystemParam.initialIntensity;%mW
num_rays=SystemParam.numLedRays; %number of total rays
num_theta=SystemParam.angleNum;%angle fibelity
radialIndex=round(num_rays/num_theta);%indexes
angularIndex=num_theta;%indexes
A=(LED_height*10^-4)*(LED_width*10^-4);
initialIntensity=P_init/A;%mW/cm2;
dA=((LED_height*10^-4)/radialIndex)*((LED_width*10^-4)/radialIndex);%area in cm^2
%origin of each ray location
Ray_Y=linspace(-LED_height/2,LED_height/2,radialIndex);%um
Ray_X=linspace(-LED_width/2,LED_width/2,radialIndex);%um
%angles of each ray at any location
alpha_ang=linspace(-ledAngle/2,ledAngle/2,angularIndex);%Rotation around X axis off Z (YZ plane rotation)
beta_ang=linspace(-ledAngle/2,ledAngle/2,angularIndex);%otation around Y axis (XZ plane rotation)
if angularIndex==1
    alpha_ang=0;
    beta_ang=0;
end
if radialIndex==1
    Ray_Y=0;
    Ray_X=0;
end

%Light intensity distribution calculations
%lambertian distribution assumtion
m=-log(2)/log(cos(ledHalfAngle));
Iray_a=cos(alpha_ang).^m;
Iray_b=cos(beta_ang).^m;
Iray_ab=Iray_a.*transpose(Iray_b);%light intensity distribution at any location

%Intensity=Iray_ab/((radialIndex^2)*sum(Iray_ab,'all'));%intensity normalized to 1 over entire LED
Intensity=ones(angularIndex,angularIndex,radialIndex,radialIndex);
%4D intensity values
for i=radialIndex
    for j=radialIndex
        Intensity(:,:,i,j)=Iray_ab;
    end
end
Intensity=Intensity/sum(sum(sum(sum(Intensity))));
% Intensity=Intensity/(sum(Intensity,'all'));%normalized 4D intensity values at x rot ang, y rot ang, x location, y location
Intensity=Intensity.*initialIntensity;%4D intensity that sums to the total initialIntensity
Int_total=sum(Intensity,'all');
end

