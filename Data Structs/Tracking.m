classdef Tracking
    properties
        %within the main model
        photon_vert_count;%counting all the times the model exits a function bc the direction is vertical
        photon_min_count;%counting all the times the model exits a function bc the direction is vertical
        transmission_count;%counting all the time the transmission case is reached
        big_dif_count_main;%counting all of the times theres a significant difference with a magnitude >10^-6
        big_dif_amt_main;%recording the significant difference amount
        big_dif_pos_main;%recording the significant difference amount
        big_dif_fun_main;%record the name of the function causing issues
        big_dif_count_meas_main;
        IO_count_main
        IO_amt_main
        IO_pos_main
        IO_fun_main
        IO_count_meas_main
            
        %within a function
        big_dif_count_func;%counting all of the times theres a significant difference with a magnitude >10^-6
        big_dif_amt_func;%recording the significant difference amount
        big_dif_pos_func;%recording the significant difference amount
        big_dif_fun_func;%record the name of the function causing issues
        big_dif_count_meas_func;
        IO_count_func;
        IO_amt_func;
        IO_pos_func;
        IO_fun_func;
        IO_count_meas_func;
        case1count;
        case2count;
        case3count;
        whiletravcount;
        
        
        
    end
end
        