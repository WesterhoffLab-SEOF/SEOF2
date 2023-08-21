function [x2,Y2_cell, num_max,Total_side] = Bins(x1,y1, incr, xlen, rad)
% Collects all data points and sums light in given increments
%emissions divided by the surface area of the fiber section to calculate per cm2
%INPUTS:: x1(:,1)= measurement location x (cm), x1(:,2)= measurement location y (cm),
%y1=intensity measured (mW), incr = fiber section size(cm), xlen = length of
%fiber (um),
SystemParam.SMA_flushlength=1*10^4;
SystemParam.SMA_totallength=2.5*10^4;
xa = x1(:,1);
xb = x1(:,2);
indexes=find(y1);
xa_short=xa(indexes,1).*10^-4;
y1_short=y1(indexes);
radius=rad*1e-4;%radius
Total_side=sum(y1,'all')
L1=length(xa);
L2=length(xa_short);
incr=incr*10^-4;

num_max = floor(((xlen/10e3)-(SystemParam.SMA_totallength/10e3))/(incr))+1;%dividing length of fiber by the increments, then adding one to have data at each end of bin
x2 = zeros(1,num_max+2);
y2 = zeros(1,num_max+2);

%first two increments are within the SMA flush length
Inc1=SystemParam.SMA_flushlength*10^-4;
Inc2=(SystemParam.SMA_totallength-SystemParam.SMA_flushlength)*10^-4;

%first two points have irregular indexes if a part of the sma connector
x2(1)=-(Inc2+Inc1);
x2(2)=-Inc1;
indX1=find(xa_short<=Inc1);
indX2_1= find(xa_short<=(Inc2+Inc1));
indX2_2=find(xa_short>(Inc1));
indX2=intersect(indX2_1,indX2_2);%logical intersection
totalY1=sum(y1_short(indX1));
totalY2=sum(y1_short(indX2));
y2(1)=totalY1/(2*pi*Inc1*radius);
y2(2)=totalY2/(2*pi*Inc2*radius);

for i=0:num_max-1
    ind2=i+3;%start at the 3rd index in the x2 and y2 vectors
    start_x2=incr*i;
    start_xa=incr*i+(-x2(1));
    stop_xa=start_xa+incr;
    x2(ind2)=start_x2;
    indxa_1=find(xa_short<=stop_xa); indxa_2=find(xa_short>start_xa);
    indxa=intersect(indxa_1,indxa_2);%logical intersection
    
    totaly2=sum(y1_short(indxa));
    y2(ind2)=totaly2/(2*pi*incr*radius);
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
x2(end)=incr*(i+1);
y2(end)=y2(end-1);
Y2_cell={y2};
end