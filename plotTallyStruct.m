function result = plotTallyStruct(Tally)

% % 
% % %%%%% this is all just extracting the data from the Tally struct to find where the biggest differences in coming from for trouble shooting
% % %%%% i didn't comment this well nor is this probably the best method to
% % %%%% extract the data. but it works for now
% % %figuring out where the gaps are, sort the strings of the difference values
% % [m1,m2,m3,m4]=size((Tally.big_dif_fun_main));%[iteration, fiber number, aa limit, number of times thing got assigned to the tally]
% % [f1,f2,f3,f4,f5,f6]=size((Tally.big_dif_fun_func));%[iteration, fiber number, aa limit,xx,yy, number of times thing got assigned to the tally]
% % 
% % maincount=numel(Tally.big_dif_count_meas_main);%total number of times we count a "big difference" in the main code
% % funcount=numel(Tally.big_dif_count_meas_func);%total number of times we count a "big difference" within each function
% % IO_maincount=numel(Tally.IO_count_meas_main);
% % IO_funcount=numel(Tally.IO_count_meas_func);
% % 
% % bigmain=strings(1,maincount);%empty vector for identifying the part of the code that has the big difference
% % dif_main=zeros(1,maincount);%empty vector for recording the big difference within the main
% % dif_fun=zeros(1,funcount);%empty vector for identifying the function that has the big difference
% % bigfun=strings(1,funcount);%empty vector for recording the big difference within the function
% % cmain=1;%index
% % cfun=1;%index for the main code
% % %for the in and out check
% % IO_bigmain=strings(1,IO_maincount);%empty vector for identifying the part of the code that has the big difference
% % IO_dif_main=zeros(1,IO_maincount);%empty vector for recording the big difference within the main
% % IO_dif_fun=zeros(1,IO_funcount);%empty vector for identifying the function that has the big difference
% % IO_bigfun=strings(1,IO_funcount);%empty vector for recording the big difference within the function
% % IO_cmain=1;%index
% % IO_cfun=1;%index for the main code
% % 
% % for i=1:m1 %iteration number
% %     for j=1:m2 %fiber number
% %         for k=1:m3 %aa number
% %             %set up emtry string vector
% %            count_inst=Tally.big_dif_count_main(i,j,k);%find maximum index where there is a location in main that had
% %            IO_count_inst=Tally.IO_count_main(i,j,k);%find maximum index where there is a location in main that had
% % 
% %            %a big difference
% %             for l=1:(count_inst)
% %                 if ischar(Tally.big_dif_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
% %                    
% %                     bigmain(1,cmain)=Tally.big_dif_fun_main{i,j,k,l}{1,1};%assign the string value
% %                     dif_main(1,cmain)=cell2mat(Tally.big_dif_count_meas_main{i,j,k,l}); %assign the difference value
% %                     cmain=cmain+1; %update the index
% %                 end
% %             end
% %             for l=1:(IO_count_inst)
% %                 if ischar(Tally.IO_fun_main{i,j,k,l}{1,1})%if we actually have something assigned here
% %                    
% %                     IO_bigmain(1,IO_cmain)=Tally.IO_fun_main{i,j,k,l}{1,1};%assign the string value
% %                     IO_dif_main(1,IO_cmain)=cell2mat(Tally.IO_count_meas_main{i,j,k,l}); %assign the difference value
% %                     IO_cmain=IO_cmain+1; %update the index
% %                 end
% %             end
% %             %assign
% %             for l=1:f4
% %                 for m=1:f5
% %                     %find each location where there is a function that had
% %                     %a big difference
% %                     count_inst=find(Tally.big_dif_fun_func{i,j,k,l,m});%
% %                     IO_count_inst=Tally.IO_count_func(i,j,k,l,m);
% % 
% %                     for o=1:length(count_inst)
% %                         n=count_inst(o);%the index
% %                             %double check the difference function isn't
% %                             %empty
% %                             if isempty(Tally.big_dif_fun_func{i,j,k,l,m,n})
% %                                 fun_str_inst="";
% %                             else
% %                                 fun_str_inst=Tally.big_dif_fun_func{i,j,k,l,m,n};
% %                                  dif_fun(1,cfun)=Tally.big_dif_count_meas_func{i,j,k,l,m,n};
% %                                  bigfun(1,cfun)=fun_str_inst;
% %                                  cfun=cfun+1;
% % 
% %                             end
% %                     end
% %                     if IO_count_inst~=0
% %                         for o=1:length(IO_count_inst)
% %                             n=IO_count_inst(o);%the index
% %                                 %double check the difference function isn't
% %                                 %empty
% %                                 if isempty(Tally.IO_fun_func{i,j,k,l,m,n})
% %                                     fun_str_inst="";
% %                                 else
% %                                     fun_str_inst=Tally.IO_fun_func{i,j,k,l,m,n};
% %                                     IO_dif_fun(1,IO_cfun)=Tally.IO_count_meas_func{i,j,k,l,m,n};
% %                                     IO_bigfun(1,IO_cfun)=fun_str_inst;
% %                                     IO_cfun=IO_cfun+1;
% % 
% %                                 end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % idfun=find(bigfun~="");
% %     if sum(dif_fun,'all')==sum(dif_fun(idfun),'all')
% %         bigfun=bigfun(idfun);
% %         fun_cats=unique(bigfun);%create a category of each unique function reported with a big dif
% %         trouble_functions=categorical(bigfun,fun_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % idmain=find(bigmain~="");
% %     if sum(dif_main,'all')==sum(dif_main(idmain),'all')
% %         bigmain=bigmain(idmain);
% %         main_cats=unique(bigmain);%create a category of each unique function reported with a big dif
% %         trouble_main=categorical(bigmain,main_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % %main_cats=unique(bigmain);%create a category of each unique function reported with a big dif
% % summary(trouble_functions)%get a summary of all of the trouble functions
% % %trouble_main=categorical(bigmain,main_cats);%main_categorical);
% % summary(trouble_main) %get a summary of all of the trouble areas in the main code
% % 
% % IO_idfun=find(IO_bigfun~="");
% %     if sum(IO_dif_fun,'all')==sum(IO_dif_fun(IO_idfun),'all')
% %         IO_bigfun=IO_bigfun(IO_idfun);
% %         IO_fun_cats=unique(IO_bigfun);%create a category of each unique function reported with a big dif
% %         IO_trouble_functions=categorical(IO_bigfun,IO_fun_cats);%,fun_categorical);
% %     else
% %         error('figure out the indexing issues')
% %     end
% % IO_idmain=find(IO_bigmain~="");%    ;
% % IO_main_cats=unique(IO_bigmain(IO_idmain));%create a category of each unique function reported with a big dif
% % summary(IO_trouble_functions);%get a summary of all of the trouble functions
% % IO_trouble_main=categorical(IO_bigmain,IO_main_cats);%main_categorical);
% % summary(IO_trouble_main) %get a summary of all of the trouble areas in the main code
% % 
% % %create a table with summary values of  issues caused by each
% % %function
% % %empty cell for all of the categories of each big_dif item
% % dif_f=cell(length(fun_cats),1);
% % dif_m=cell(length(main_cats),1);
% % IO_dif_f=cell(length(fun_cats),1);
% % IO_dif_m=cell(length(main_cats),1);
% % %empty vectory for all of the categories of each big_dif item
% % dif_contrib_m=zeros(length(main_cats),2);
% % dif_contrib_f=zeros(length(fun_cats),2);
% % 
% % IO_dif_contrib_m=zeros(length(IO_main_cats),2);
% % IO_dif_contrib_f=zeros(length(IO_fun_cats),2);
% % for i=1:length(main_cats)%for each categorical in the main function
% %     dif_m(1,i)={dif_main(1,(trouble_main==main_cats(i)))};% max a cell 
% %     dif_contrib_m(i,1)=sum(cell2mat(dif_m(1,i)),'all');
% %     dif_contrib_m(i,2)=sum(abs(cell2mat(dif_m(1,i))),'all');
% % end
% % for i=1:length(fun_cats)
% %     dif_f(1,i)={dif_fun(1,(trouble_functions==fun_cats(i)))};
% %     dif_contrib_f(i,1)=sum(cell2mat(dif_f(1,i)),'all');
% %     dif_contrib_f(i,2)=sum(abs(cell2mat(dif_f(1,i))),'all');
% % end
% % 
% % for i=1:length(IO_main_cats)%for each categorical in the main function
% %     IO_dif_m(1,i)={IO_dif_main(1,(IO_trouble_main==IO_main_cats(i)))};% max a cell 
% %     IO_dif_contrib_m(i,1)=sum(cell2mat(IO_dif_m(1,i)),'all');
% %     IO_dif_contrib_m(i,2)=sum(abs(cell2mat(IO_dif_m(1,i))),'all');
% % end
% % for i=1:length(IO_fun_cats)
% %     IO_dif_f(1,i)={IO_dif_fun(1,(IO_trouble_functions==IO_fun_cats(i)))};
% %     IO_dif_contrib_f(i,1)=sum(cell2mat(IO_dif_f(1,i)),'all');
% %     IO_dif_contrib_f(i,2)=sum(abs(cell2mat(IO_dif_f(1,i))),'all');
% % end
% % 
% % 
% % Changed=dif_contrib_m(:,1);%difference value vector for the table
% % Total=dif_contrib_m(:,2);%absolute value vector for the table
% % Function=transpose(main_cats); %function
% % TabM=table(Changed,Total,'RowNames',Function)
% % Changed=dif_contrib_f(:,1);%difference value vector for the table
% % Total=dif_contrib_f(:,2);%absolute value vector for the table
% % Function=transpose(fun_cats);
% % TabF=table(Changed,Total,'RowNames',Function)
% % 
% % 
% % IO_Changed=IO_dif_contrib_m(:,1);%difference value vector for the table
% % IO_Total=IO_dif_contrib_m(:,2);%absolute value vector for the table
% % IO_Function=transpose(IO_main_cats); %function
% % IO_TabM=table(IO_Changed,IO_Total,'RowNames',IO_Function)
% % IO_Changed=IO_dif_contrib_f(:,1);%difference value vector for the table
% % IO_Total=IO_dif_contrib_f(:,2);%absolute value vector for the table
% % IO_Function=transpose(IO_fun_cats);
% % IO_TabF=table(IO_Changed,IO_Total,'RowNames',IO_Function)
% % 

end