function Tally = initializeTally(it_num, SystemParam, numRandomIters, a, b)
    %empty storage for a numeric count of some stuff for trouble
    %shooting
    Tally.photon_vert_count=zeros(it_num,SystemParam.numFibers,numRandomIters);
    Tally.minPhotons_count=zeros(it_num,SystemParam.numFibers,numRandomIters);
    Tally.transmission_count=zeros(it_num,SystemParam.numFibers,numRandomIters);%all the time the transmission count is reached
    Tally.big_dif_count_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%counting all of the times theres a significant difference with a magnitude >10^-6
    Tally.big_dif_amt_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%recording the significant difference amount
    Tally.big_dif_pos_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%recording the significant difference amount
    Tally.big_dif_fun_main=cell(it_num,SystemParam.numFibers,numRandomIters);%record the name of the function causing issues
    Tally.big_dif_count_meas_main=cell(it_num,SystemParam.numFibers,numRandomIters);
    Tally.IO_count_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%counting all of the times theres a significant difference with a magnitude >10^-6
    Tally.IO_amt_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%recording the significant difference amount
    Tally.IO_pos_main=zeros(it_num,SystemParam.numFibers,numRandomIters);%recording the significant difference amount
    Tally.IO_fun_main=cell(it_num,SystemParam.numFibers,numRandomIters);%record the name of the function causing issues
    Tally.IO_count_meas_main=cell(it_num,SystemParam.numFibers,numRandomIters);
    
    %within a function
    Tally.big_dif_count_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%counting all of the times theres a significant difference with a magnitude >10^-6
    Tally.big_dif_amt_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%recording the significant difference amount
    Tally.big_dif_pos_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%recording the significant difference amount
    Tally.big_dif_fun_func=cell(it_num,SystemParam.numFibers,numRandomIters,a,b);%record the name of the function causing issues
    Tally.big_dif_count_meas_func=cell(it_num,SystemParam.numFibers,numRandomIters,a,b);
    Tally.IO_count_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%counting all of the times theres a significant difference with a magnitude >10^-6
    Tally.IO_amt_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%recording the significant difference amount
    Tally.IO_pos_func=zeros(it_num,SystemParam.numFibers,numRandomIters,a,b);%recording the significant difference amount
    Tally.IO_fun_func=cell(it_num,SystemParam.numFibers,numRandomIters,a,b);%record the name of the function causing issues
    Tally.IO_count_meas_func=cell(it_num,SystemParam.numFibers,numRandomIters,a,b);
    Tally.case1count=zeros(it_num,SystemParam.numFibers,numRandomIters);
    Tally.case2count=zeros(it_num,SystemParam.numFibers,numRandomIters);
    Tally.case3count=zeros(it_num,SystemParam.numFibers,numRandomIters);
    Tally.whiletravcount=zeros(it_num,SystemParam.numFibers,numRandomIters);
end