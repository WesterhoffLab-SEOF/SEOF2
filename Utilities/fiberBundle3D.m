function [circles,combocircles] = fiberBundle3D(num_fib,fiber_D)
%%%%%%%%takes the number of fibers, the fiber diameter of each (assumed the
%%%%%%%%same), then reports the center of each fiber [circles] and a black
%%%%%%%%+white image showing the packing geometry. Uses the minimum circle
%%%%%%%%packing diameter (NS note: don't remember how I did this. it works
%%%%%%%%consistently so am not trying to re-learn to edit.)
if num_fib==1
           %number of fibers
            d_fib=fiber_D;%(um)
            r_fib=d_fib/2;%(um)
            A_fib=num_fib*pi*r_fib^2;%(um^2)
            
            %SMA
            %optimized smallest diameter circle for n unit circles to fit into..
            d_SMA=(d_fib)*(1);%(um) https://mathworld.wolfram.com/CirclePacking.html
            r_SMA=d_SMA/2;%(um)
            A_SMA=pi*r_SMA^2;%(um^2)
            
            %create 2D logical Image of circles
            imageSizeX=2500;%
            imageSizeY=2500;%
            %each pixel is 1x1um
            [column_image,row_image]=meshgrid(1:imageSizeX,1:imageSizeY);
            center_image=[round(imageSizeX/2),round(imageSizeY/2)];
            circles=(row_image-(center_image(2))).^2+(column_image-(center_image(1))).^2 <=r_fib.^2;
            combocircles=circles;

        elseif num_fib==4
            num_fib=4;%number of fibers
            d_fib=fiber_D;%(um)
            r_fib=d_fib/2;%(um)
            A_fib=num_fib*pi*r_fib^2;%(um^2)
            
            %SMA
            %optimized smallest diameter circle for n unit circles to fit into..
            d_SMA=(d_fib)*(1+sqrt(2));%(um) https://mathworld.wolfram.com/CirclePacking.html
            r_SMA=d_SMA/2;%(um)
            A_SMA=pi*r_SMA^2;%(um^2)
            
            %Orientation
            %choose
            %or_fib1=0;%fiber centered on 0 radians
            or_fib1=pi/4;%fiber centered on pi/4 radians
            
            %Fiberlocationfromcenter
            a_fibs=[or_fib1,or_fib1+(pi/2),or_fib1+(pi),or_fib1+(3*pi/2)];%angles from x-axis
            cent_fib=((d_SMA-(2*d_fib))/2)+r_fib;%angles from y axis
            %XY coordinates
            XY1=round(cent_fib.*[cos(a_fibs(1)),sin(a_fibs(1))]);
            XY2=round(cent_fib.*[cos(a_fibs(2)),sin(a_fibs(2))]);
            XY3=round(cent_fib.*[cos(a_fibs(3)),sin(a_fibs(3))]);
            XY4=round(cent_fib.*[cos(a_fibs(4)),sin(a_fibs(4))]);

            %create 2D logical Image of circles
            imageSizeX=2500;%
            imageSizeY=2500;%
            %each pixel is 1x1um
            [column_image,row_image]=meshgrid(1:imageSizeX,1:imageSizeY);
            center_image=[round(imageSizeX/2),round(imageSizeY/2)];
            
            circles=zeros(imageSizeX,imageSizeY,4);
            circles(:,:,1)=(row_image-(center_image(2)+XY1(2))).^2+(column_image-(center_image(1)+XY1(1))).^2 <=r_fib.^2;
            circles(:,:,2)=(row_image-(center_image(2)+XY2(2))).^2+(column_image-(center_image(1)+XY2(1))).^2 <=r_fib.^2;
            circles(:,:,3)=(row_image-(center_image(2)+XY3(2))).^2+(column_image-(center_image(1)+XY3(1))).^2 <=r_fib.^2;
            circles(:,:,4)=(row_image-(center_image(2)+XY4(2))).^2+(column_image-(center_image(1)+XY4(1))).^2 <=r_fib.^2;
            combo12=or(circles(:,:,1),circles(:,:,2));
            combo34=or(circles(:,:,3),circles(:,:,4));
            combocircles=or(combo12,combo34);

        elseif num_fib==19
            num_fib=19;%number of fibers
            d_fib=fiber_D;%(um)
            r_fib=d_fib/2;%(um)
            A_fib=num_fib*pi*r_fib^2;%(um^2)
            
            %SMA
            %optimized smallest diameter circle for n unit circles to fit into..
            r_SMA=r_fib*(1+sqrt(2)+sqrt(6));%(um) https://mathworld.wolfram.com/CirclePacking.html
            d_SMA=r_fib*2;%(https://link.springer.com/content/pdf/10.1023/A:1005091317243.pdf
            A_SMA=pi*r_SMA^2;%(um^2)
            
            %create 2D logical Image of circles
            imageSizeX=2500;%
            imageSizeY=2500;%
            %each pixel is 1x1um
            [column_image,row_image]=meshgrid(1:imageSizeX,1:imageSizeY);
            center_image=[round(imageSizeX/2),round(imageSizeY/2)];
       
            
            circles=zeros(imageSizeX,imageSizeY,19);
            %center circle
            circles(:,:,1)=(row_image-(center_image(2))).^2+(column_image-(center_image(1))).^2 <=r_fib.^2;
            combocircles=circles(:,:,1);
            %cicles 2-7 immediately surrounding the center circle
            center_fib27=d_fib;
            for i=1:6
                combocirci=combocircles;%initialize the combination vector
                a_fibs=(i-1)*(2*pi)/6;%angle of fiber from the center
                XY=round(center_fib27.*[cos(a_fibs),sin(a_fibs)]);
                circles(:,:,i+1)=(row_image-(center_image(2)+XY(2))).^2+(column_image-(center_image(1)+XY(1))).^2 <=r_fib.^2;
                combocircles=or(combocirci,circles(:,:,i+1));
            end
            %for fibers 8thru 19
             center_fib819=sqrt((2*(d_fib^2))-(2*d_fib*d_fib*cos(pi-(2*pi/24))));%solution done w law of cosines
            for i=1:12
                combocirci=combocircles;%initializw the combo vecotr
                a_fibs=((i-1)*(2*pi)/12)+(pi/12);%angle of fiber from the center
                XY=round(center_fib819.*[cos(a_fibs),sin(a_fibs)]);
                circles(:,:,i+7)=(row_image-(center_image(2)+XY(2))).^2+(column_image-(center_image(1)+XY(1))).^2 <=r_fib.^2;
                combocircles=or(combocirci,circles(:,:,i+7));
            end


end

end

