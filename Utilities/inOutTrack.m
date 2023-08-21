%trouble shooting function 
function [Tally,Differenceamount,Diffamountpos] =  inOutTrack(In_Pow,Out_Pow,SystemParam,name_fun,main_or_func,Global_Index,Tally)%y,gg)
    %simple function that takes a scalar of the input power, and a vector
    %of the various output powers,the name of the function, then find is the difference is
    %significantly large. Name must be a string.
    %Global_Index=[iteration,h,aa,xx,yy];
    %main_or_func 0=main function, 1 = within a function called
    st_pow=sum(In_Pow,'all');
    out_pow=sum(Out_Pow,'all');
    iteration=Global_Index(1);h=Global_Index(2);aa=Global_Index(3);xx=Global_Index(4);yy=Global_Index(5);
    Differenceamount = sum(In_Pow,'all') - sum(Out_Pow,'all');
    Diffamountpos = abs(Differenceamount);
    %include a IO_amt_counted with counted cells
    if main_or_func==0%in the main function
        if st_pow<out_pow
            Tally.IO_count_main(iteration,h,aa)=Tally.IO_count_main(iteration,h,aa)+1;
            Tally.IO_amt_main(iteration,h,aa)=Tally.IO_amt_main(iteration,h,aa)+Differenceamount;
            Tally.IO_pos_main(iteration,h,aa)=Tally.IO_pos_main(iteration,h,aa)+Diffamountpos;
            Tally.IO_fun_main{h,iteration,aa,Tally.IO_count_main(iteration,h,aa)}={name_fun};
            Tally.IO_count_meas_main{h,iteration,aa,Tally.IO_count_main(iteration,h,aa)}={Differenceamount};
        end
    elseif main_or_func==1%within a function
        if st_pow<out_pow
            Tally.IO_count_func(iteration,h,aa,xx,yy)=Tally.IO_count_func(iteration,h,aa,xx,yy)+1;
            Tally.IO_amt_func(iteration,h,aa,xx,yy)=Tally.IO_amt_func(iteration,h,aa,xx,yy)+Differenceamount;
            Tally.IO_pos_func(iteration,h,aa,xx,yy)=Tally.IO_pos_func(iteration,h,aa,xx,yy)+Diffamountpos;
            Tally.IO_fun_func(iteration,h,aa,xx,yy,Tally.IO_count_func(iteration,h,aa,xx,yy))={name_fun};
            Tally.IO_count_meas_func(iteration,h,aa,xx,yy,Tally.IO_count_func(iteration,h,aa,xx,yy))={Differenceamount};
        end
    end
end