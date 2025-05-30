function [x2,Y2_cell, num_max,Total_side] = Bins(x1,y1, incr, xlen, rad,SystemParam)
% Collects all data points and sums light in given increments
%emissions divided by the surface area of the fiber section to calculate per cm2
%INPUTS:: x1(:,1)= measurement location x (cm), x1(:,2)= measurement location y (cm),
%y1=intensity measured (mW), incr = fiber section size(cm), xlen = length of
%fiber (um),

xa = x1(:,1);
indexes=find(y1);
xa_short=xa(indexes,1).*10^-4;%x position in cm
y1_short=y1(indexes);
incr=incr*10^-4;%increment in cm
radius=rad*1e-4;%radius in cm

%if there's an SMA connector, the x vector will record differently (x
%within the sma connect is negative, the 0 starts outside the SMA
%connector)
if SystemParam.SMA==1
    num_max = floor(((xlen/10e3)-(SystemParam.smaTotalLength/10e3))/(incr))+1;%dividing length of fiber by the increments, then adding one to have data at each end of bin
    x2 = zeros(1,num_max+2);
    y2 = zeros(1,num_max+2);
    
    %first two measurement points will be within the SMA flush length
    Inc1=SystemParam.smaFlushLength*10^-4;%in cm
    Inc2=(SystemParam.smaTotalLength-SystemParam.smaFlushLength)*10^-4;%in cm
    
    %first two points have irregular indexes if a part of the sma connector
    x2(1)=-(Inc2+Inc1);
    x2(2)=-Inc1;
    indX1=find(xa_short<=Inc1);
    indX2_1= find(xa_short<=(Inc2+Inc1));
    indX2_2=find(xa_short>(Inc1));
    indX2=intersect(indX2_1,indX2_2);%logical intersection
    totalY1=sum(y1_short(indX1));
    totalY2=sum(y1_short(indX2));
    y2(1)=totalY1;%/(2*pi*Inc1*radius);
    y2(2)=totalY2;%/(2*pi*Inc2*radius);
    start_pos=3;%the initial vector;
    adj_for_sma=1;%don't adjust the start and stop indexes for the sma connector
    %xa_short=xa_short-x2(1);%make everything -2.5cm
else%if no sma connector
    num_max = floor((xlen/10e3)/(incr))+1;%dividing length of fiber by the increments, then adding one to have data at each end of bin
    x2 = zeros(1,num_max);
    y2 = zeros(1,num_max);
    start_pos=1;
    adj_for_sma=0;%don't adjust the start and stop indexes for the sma connector
end


Total_side=sum(y1,'all');%.*2*pi*radius*xlen/10e3;%power in W


%setting the x at the "0 measurement)
x2(start_pos)=0+(incr/2);
start_xa=x2(start_pos)-(incr/2)+(adj_for_sma*(Inc2+Inc1));
stop_xa=start_xa+incr;
indxa_1=find(xa_short<=stop_xa); indxa_2=find(xa_short>start_xa);
indxs=intersect(indxa_1,indxa_2);%logical intersection
totaly2=sum(y1_short(indxs));
y2(start_pos)=totaly2;%/(2*pi*incr*radius);

for i=1:num_max-1
    center_x2=incr*i;%center on the whole number
    start_xa=incr*i-(incr/2)+(adj_for_sma*(Inc2+Inc1));
    stop_xa=start_xa+incr;
    
    
    ind2=i+start_pos;%select index. start at the 3rd index in the x2 and y2 vectors
    x2(ind2)=center_x2;% x vec is in the center of the measurement
    indxa_1=find(xa_short<=stop_xa); indxa_2=find(xa_short>start_xa);
    indxa=intersect(indxa_1,indxa_2);%logical intersection
    
    totaly2=sum(y1_short(indxa));
    y2(ind2)=totaly2;%/(2*pi*incr*radius);
end


% for i=0:num_max-1
%     start = incr*i;
%     stop = start+incr;
%     total = 0;
%     for a=1:length(xa_short)
%         if (xa_short(a)>=start) && (xa_short(a) <= stop) && (y1_short(a) > 0)
%             total = total + y1_short(a);
%         end
%     end
%     x2(i+1) = start;
%     y2(i+1) = total/(2*pi*incr*radius);% sum of light divided by surface area of section
% end
%setting final value

x2(end)=xlen/10e3-(incr/2)-(adj_for_sma*(Inc2+Inc1));%want value 0.5cm from the tip of the fiber, if we adjust for sma all x values sub -2.5cm
start_xa=x2(end)-(incr/2);
stop_xa=start_xa+incr;
indxa_1=find(xa_short<=stop_xa); indxa_2=find(xa_short>start_xa);
indxa=intersect(indxa_1,indxa_2);%logical intersection
totaly2=sum(y1_short(indxa));
y2(end)=totaly2;%/(2*pi*incr*radius);

Y2_cell={y2};
end