function[theta_i,theta_t,theta_c,theta_ih,theta_rh,V_transmit,V_reflect,direction] = Snells(Vin,nhat,ni1,nt1,horz_surf,direction)
%calculates the snell's law different angles, gives back a new light
%direction%vectors
ni=real(ni1);
nt=real(nt1);
V=Vin;
if V(2)==0
    dirmult=1;
else
    dirmult=sign(V(2));
end
if length(nhat)==3
    nhat=nhat(1:2);%z direction is never used in the 2d model but just in case
end
if(horz_surf==1)%impacting horizontal surface
    nv=[0,dirmult*-1];%normal vector is opposite the y dir of the impacting 
    Rcw90=[0,1;-1,0];%90 degree clockwiserotation matrix
    n_hor=transpose(Rcw90*transpose(nv));%the horizontal vector is 90 deg rotated from the vertical vector
    n_=nhat.*[1,dirmult];%make sure it's going opposite of dir vector (convention is to set it as [dy/dx,-1], or [-1,dx/dy]
    rotcalc=asin((cross([nv,0],[n_,0]))/(norm(n_)*norm(nv)));
    theta_rot=rotcalc(3);
else
    n_hor=[-1*sign(Vin(1)),0];
    Rccw90=[0,-1;1,0];
    nv=transpose(Rccw90*transpose(n_hor));%the vertical vector is 90 deg rotated from the horizontal vector
    n_=nhat.*[sign(Vin(1)),1];%make sure it's going opposite of dir vector (convention is to set it as [dy/dx,-1], or [-1,dx/dy]
    rotcalc=asin((cross([n_hor,0],[n_,0]))/(norm(n_)*norm(n_hor)));
    theta_rot=rotcalc(3);
end
%rotation matrices
Rccw=[cos(theta_rot) ,-sin(theta_rot);sin(theta_rot), cos(theta_rot)];%counter clock wise (in correct dir)
Rcw=[cos(theta_rot) ,sin(theta_rot);-sin(theta_rot), cos(theta_rot)];%clock wise (opposite dir to put things back in place)
%n_hor=[sign(Vin(1))*-1,0,0]%horizontal vector from which the Velocities are calculated
if length(Vin)==2 %make sure the vectors are appropriate size for a cross prouct
    V=[Vin,0];
end

theta_ih=abs(asin(norm(cross(V,[n_hor,0])/(norm(V)*norm([n_hor,0])))))*sign(V(1))*sign(V(2));%will assign direction based on Vin


theta_i=abs(asin(norm(cross(V,[n_,0])/(norm(V)*norm([n_,0])))));%will assign direction based on Vin

%%snells law%%
theta_c=asin(nt/ni);
theta_t = asin((ni/nt)*sin(theta_i));

%rotate the incoming velocity vector
Vinp=transpose(Rccw*transpose(V(1:2)));%rotated incoming vector
if Vinp(2)==0
    dirmultp=1;
else
    dirmultp=sign(V(2));
end
if(horz_surf==1)%impacting horizontal surface
    Vrefp=abs(Vinp).*[sign(Vinp(1)),-1*dirmultp];%changes dir in the y' direction
    if ((theta_i>theta_c)&& isreal(theta_c)) || ~isreal(theta_t)
        Vtransp=[sign(Vinp(1))*1,0];
    else
        Vtransp=[sign(Vinp(1))*cos(pi/2-theta_t),dirmultp*sin(pi/2-theta_t)];%need to get thetat from x' coord
    end
else%impacting vertical surface
    Vrefp=abs(Vinp).*[-1*sign(Vinp(1)),dirmultp];%changes dir in the x' direction

    if ((theta_i>theta_c)&& isreal(theta_c)) || ~isreal(theta_t)
        Vtransp=[0,dirmultp*1];
    else
        Vtransp=[sign(Vinp(1))*cos(theta_t),dirmultp*sin(theta_t)];%theta t from x' coordinate already
    end
end
%Rotate back to normal
    %V_in_check=transpose(Rcw*transpose(Vinp));
    V_reflect=transpose(Rcw*transpose(Vrefp));
   
    V_transmit=transpose(Rcw*transpose(Vtransp));
    if sign(V_reflect(1))*sign(Vin(1))==-1
        direction=-1*direction;
    end
    %angle of the V_reflect vector off of the horizontal
    theta_rh=sign(V_reflect(2))*abs(asin(norm(cross([V_reflect,0],[direction,0,0])/(norm([V_reflect,0])*norm([direction,0,0])))));%+/- based on V_reflect ydir
    
end
