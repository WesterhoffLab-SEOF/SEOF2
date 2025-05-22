function [I_scatter,theta_scatter]=scatter_cone(SystemParam,theta_main,vert_surf)
%create a vector for the scatter cone of the relative magnitude of each scattered
%intensity, and the angles off of the x axis, and the velocity of each ray
%this function ONLY for vertical surfaces
if vert_surf==0
    error('the scatter cone function is only for vertical surfaces')
end

scatterNum=SystemParam.scatterNum;%maximum number of scatter rays to follow
max_angle=SystemParam.maxScatterAngle;%maximum angle off of the original to create the scatter cone
lim_angle=SystemParam.ledAngle/2;%limiting angle

if scatterNum==1 || scatterNum==0 || abs(theta_main)>=lim_angle%if theres only 1 or no scatter angles,oor the angle is rly sharp, everything stays the same
    theta_scatter=theta_main;
    %V_scatter=V_main;
    I_scatter=1;
else%if there are multiple scattering rays and the input angle is reasonable to scatter from
    theta_scatter=linspace(theta_main-max_angle,theta_main+max_angle,scatterNum);%create a vector with the max
    %check that no scattering angle is more than the limit (i.e. would be
    %going outside of the fiber)
    lessthanmin=find(theta_scatter<-lim_angle);%index of the less than min
    morethanmax=find(theta_scatter>lim_angle);%index of the more than max
    if ~isempty(lessthanmin)%if there are scattering angles less than the minimm
        if length(lessthanmin)==1%if just the first angle in the vector (min angle)
            theta_scatter(1)=-lim_angle;%replace that abngle with the limit;
            % %consider creating a new scatter cone between -lim_angle and
            %the theta_main+dif_lim
            %dif_lim=abs(-lim_angle-theta_main);%find the difference between the limit and the main angle
            %theta_scatter=linspace(-lim_angle,theta_main+dif_lim,scatterNum);
        else
            index_threshold=max(lessthanmin);%find the maximum index (ie the one least furthest off of the limit
            theta_scatter(index_threshold)=-lim_angle;%replace that one with the limit
            length_ones=index_threshold-1;
            theta_scatter(1:index_threshold-1)=100*ones(1,length_ones);%everything that was below the limit is replaced with 100, which will be removed
            %            %consider creating a new scatter cone between -lim_angle and
            %the theta_main+dif_lim
            %dif_lim=abs(-lim_angle-theta_main);%find the difference between the limit and the main angle
            %theta_scatter=linspace(-lim_angle,theta_main+dif_lim,scatterNum);
        end
    end
    if ~isempty(morethanmax)
        if length(morethanmax)==1
            theta_scatter(morethanmax)=lim_angle;
            %            %consider creating a new scatter cone between theta_main-lim_angle and
            %the theta_main+dif_lim
            %dif_lim=abs(lim_angle-theta_main);%find the difference between the limit and the main angle
            %theta_scatter=linspace(theta_main-dif_lim,lim_angle,scatterNum);
        else%if length(morethanmax)>1%if there are multiple more than max
            index_threshold=min(morethanmax);%find the index closest to the limit
            length_ones=length(theta_scatter)-index_threshold;
            theta_scatter(index_threshold)=lim_angle;


        theta_scatter(1,index_threshold+1:end)=100*ones(1,length_ones);%everything that above below the limit is replaced with 100, which will be removed
        %            %consider creating a new scatter cone between -lim_angle and
        %the theta_main+dif_lim
        %dif_lim=abs(-lim_angle-theta_main);%find the difference between the limit and the main angle
        %theta_scatter=linspace(-lim_angle,theta_main+dif_lim,scatterNum);
        end
    end
    theta_scatter=theta_scatter(theta_scatter~=100);%select for just the non 100 values
    
    theta_cone=theta_scatter-theta_main;%angle off of the initial theta for use in intensity calculations
    I_temp=1+cos(theta_cone).^2;%lambertian distribution of intensity
    sum_Itemp=sum(I_temp);%sum the whole value
    I_scatter=I_temp./sum_Itemp;%normalize to the total
    %V_scatter=[sign(V_main(1))*abs(cos(theta_scatter)),sin(theta_scatter)];%going to be going in the same x direction as the initial/main ray
    
end
end
