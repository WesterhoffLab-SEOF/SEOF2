function [Ent_Int,Theta,y0,perc_hit,perc_ent,Entering_Int,Entering_angle,Entering_X,Entering_Y,incoming_int,incoming_ang] = enterLight3D(SystemParam,r_fib,led_d,alpha_ang,beta_ang,Ray_X,Ray_Y,Intensity_mat,circles)
%variables from system parameters
nt = SystemParam.n1;              % RI of Quartz Optical Fiber %refracted index
ni = SystemParam.n2;              % RI of "gap"/air %air index
%I_init = SystemParam.I_init;          % LED intensity, 50 mW
        nhat=[-1 0];%surface of fiber is the flat cut end
        horz_surf=0;%not a horizontal surface
        direction=1;
        %Ray_X and Ray_Y use (um), as does circles
%r_fiber=SystemParam.rfiber; 
k_air = SystemParam.kair*10^-4; %estimated light attenuation constant (1/cm)->(1/um) https://thesis.library.caltech.edu/3249/1/Baum_wa_1950.pdf
LED_dist=led_d;%*10^-3;%;mm
[a,b,c,d]=size(Intensity_mat);
[e,f,g]=size(circles);
r_fiber=r_fib;%*10^-3;%um>mm
imageSizeX=e;
imageSizeY=f;
center_image=[round(imageSizeX/2),round(imageSizeY/2)];
I_init=sum(Intensity_mat,'all');
%empty vector set ups
Ent_Int=zeros(a,c,g);
Theta=zeros(a,c,g);
y0=zeros(a,c,g);
light_entering=zeros(1,g);
incoming_int=zeros(a,c,g);
incoming_ang=zeros(a,c,g);
perc_ent=zeros(1,g);
perc_hit=zeros(1,g);



%for each fiber
for q=1:g
%empty vector set ups for within each fiber
    Entering_Int=zeros(a,b,c,d);
    Inc_Int=zeros(a,b,c,d);
    Loss_Inc=zeros(a,b,c,d);
    Entering_angle=zeros(a,b,c,d);
    Entering_X=zeros(a,b,c,d);
    Entering_Y=zeros(a,b,c,d);
    Entering_rad=zeros(a,b,c,d);
    %calculate launched into fiber values
    for i=1:c%Ray X
        for j=1:d%Ray Y
            for k=1:a%angle alpha
                for l=1:b%angle beta
                    Pos_i_ray=[Ray_X(i),Ray_Y(j),0];%initial position vector(um);
                    %transformation matrixes
                    rot_x=[1,0,0;0,cos(alpha_ang(k)),-sin(alpha_ang(k));0,sin(alpha_ang(k)),cos(alpha_ang(k))];
                    rot_y=[cos(beta_ang(l)),0,sin(beta_ang(l));0,1,0;-sin(beta_ang(l)),0,cos(beta_ang(l))];
                    rot_z=eye(3);%no rotation in Z, turns into eye matrix
                    R=rot_z*rot_y*rot_x;%ZYX tait brian rotation convention
                    dir_vec=R*[0;0;-1];%direction vector, rotated in both x and then y
                    %new position upon entering calculation
                    Znew=LED_dist;
                    Ynew=round(dir_vec(2)*(Znew/dir_vec(3)))+Ray_Y(j);
                    Xnew=round(dir_vec(1)*(Znew/dir_vec(3)))+Ray_X(i);
                    %Inc_Int(k,l,i,j)=%Intensity_mat(k,l,i,j)./(4*pi*norm([Xnew,Ynew,Znew]-Pos_i_ray).^2);
                    Loss_Inc(k,l,i,j)=exp(-k_air*norm([Xnew,Ynew,Znew]-Pos_i_ray));
                    Inc_Int(k,l,i,j)=Intensity_mat(k,l,i,j).*Loss_Inc(k,l,i,j);%norm([Xnew,Ynew,Znew]-Pos_i_ray)*10^-1);%light intensity after attenuation in air
                    
                    %check if the light hits a fiber by using the boolean
                    %image
                    %transform x and y coordinates of light vector to boolean
                    %image
                    X_trans=Xnew+center_image(1);
                    Y_trans=Ynew+center_image(2);
                    if round(X_trans)==0 %image pixel index starts at 1
                        X_trans=1;
                    elseif round(Y_trans)==0 %image pixel index starts at 1
                        Y_trans=1;
                    end
                    if (X_trans>0 && X_trans<=e) && (Y_trans>0 && Y_trans<=f) %if both of the position values are within the image range (outside also things wont hit_
                    if circles(round(X_trans),round(Y_trans),q)==1
                        %3d storage
                        light_entering(1,q)=light_entering(1,q)+Inc_Int(k,l,i,j);
                        Entering_Int(k,l,i,j)=Inc_Int(k,l,i,j);
                        Entering_angle(k,l,i,j)=acos(dot(dir_vec,transpose([0,0,1])));%entering angle off of normal
                        Entering_X(k,l,i,j)=Xnew;
                        Entering_Y(k,l,i,j)=Ynew;
                         %flattening to 2-D
                        Entering_rad(k,l,i,j)=sqrt(Xnew^2+Ynew^2);
                        
                    end
                    end
                end
            end
        end
    end
    total_entering_int=sum(Entering_Int,'all');
    light_ent_check=light_entering(1,q);
    %flattening to 2D for each fiber
    for i=1:c
        for k=1:a
         incoming_int(k,i,q)=sum(sum((Entering_Int(k,:,i,:))));%total amount of intensity entering
         incoming_ang(k,i,q)=sum(sum((Entering_angle(k,:,i,:))))/(b*d);%average angle
%         incoming_int(k,i,q)=sum(Entering_Int(k,:,i,:),'all');%total amount of intensity entering
%         incoming_ang(k,i,q)=sum(Entering_angle(k,:,i,:),'all')/(b*d);%average angle
        vi=[cos(incoming_ang(k,i,q)),sin(incoming_ang(k,i,q))];%entering direction vector
        [theta_i,theta_t,theta_c,theta_ih,~,~,~,~] = Snells(vi,nhat,ni,nt,horz_surf,direction);
        [~,Ent_Int(k,i,q)]=FresnelEq(incoming_int(k,i,q),SystemParam,theta_i,theta_t,theta_c,theta_ih,ni,nt,horz_surf);
        y0(k,i,q)=sign(Ray_X(i))*sum(sum(Entering_rad(k,:,i,:)))/(b*d);%average y0
        Theta(k,i,q)=sign(alpha_ang(k))* theta_t;
        end
    end
    total_act_Ent_Int=sum(Ent_Int(:,:,q),'all');
    I_init_check=I_init;
     perc_hit(q)=total_entering_int/I_init;%perc_hit(q)=sum(sum(incoming_int(:,:,q)))/I_init;
     perc_ent(q)=total_act_Ent_Int/I_init;%sum(Ent_Int(:,:,q),'all')/I_init;  %%    perc_ent(q)=sum(sum(Ent_Int(:,:,q)))/I_init;    
end


    


                    
                    
end

