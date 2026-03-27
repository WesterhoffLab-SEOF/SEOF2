function res = writeResults(SystemParam,iterParamfields, description, filename, sheetNum, iteration, it_num,h, FibIt, c,topoffile,randvar,itdif,xlen2true)
%write results of simulation to organized excel table that will be easier
%to import for data analysis

%current organization of data is deisgned for variable analysis for Shapiro et. al, 2025 paper 
% by running 12.5 (air and water) and 50.5cm (air only) fiber simulations for each variable set


for jj=1:(randvar+1)
    %if randvar=0 then this wont change anything
if c>1 %if multiple fibers
    if it_num>1%if there are multiple iterations
        %put all the fibers within a single iteration in the same sheet,
        %move the iterations to different sheets
        if randvar==0
            sheetnum=iteration;
        else
            sheetnum=iteration*(randvar+1)-1+(jj-1);
        end
    linenum=7+(h-1);
    if sheetnum>1
        writecell(topoffile,filename,'Sheet',sheetnum,'Range',"A1");
    end
    else
        linenum=7+(h-1);
        sheetnum=sheetNum;
    end
    writematrix('fiber no.',filename,'Sheet',sheetnum,'Range',"B6");
    writematrix(h,filename,'Sheet',sheetnum,'Range',"B"+num2str(linenum));
    description=description + " fiber no. " +h;
else
    sheetnum=sheetNum+(jj-1);
    linenum=7+(iteration-1);
    writematrix(iteration,filename,'Sheet',sheetnum,'Range',"B"+num2str(linenum));
end
Paramcurrent=cell(1,length(iterParamfields));   
for v=1:length(iterParamfields)
Paramcurrent{v}=SystemParam.(iterParamfields{v});
end
    writecell({filename},filename,'Sheet',sheetnum,'Range',"A"+num2str(linenum))
    writecell(Paramcurrent,filename,'Sheet',sheetnum,'Range',"C"+num2str(linenum));
    writematrix(imag(SystemParam.n1),filename,'Sheet',sheetnum,'Range',"F"+num2str(linenum))
    %summary data
    writematrix(FibIt(jj).pow_side_waterst(iteration,h),filename,'Sheet',sheetnum,'Range',"U"+num2str(linenum))
    writematrix(FibIt(jj).transmitted(iteration,h),filename,'Sheet',sheetnum,'Range',"V"+num2str(linenum))
    if SystemParam.waterInterface==1
        newlinenum=linenum-itdif(2);
        if newlinenum<1
            newlinenum=linenum;
        end
        writematrix(FibIt(jj).pow_side_waterst(iteration,h),filename,'Sheet',sheetnum,'Range',"W"+num2str(newlinenum))
    end

    %write the Ee(x) values to the files
    if SystemParam.SMA==1
            EeX0end0=cell2mat(FibIt(jj).Y(iteration,h));
            EeX0end=EeX0end0((floor(SystemParam.smaTotalLength*10^-4)+1):end);
    elseif SystemParam.SMA==0
            EeX0end=cell2mat(FibIt(jj).Y(iteration,h));

    end
    if SystemParam.xLen==50.5*10^4
        %write to the 50.5 start location (AM6)
        writematrix(EeX0end,filename,'Sheet',sheetnum,'Range',"AM"+num2str(linenum))
        if xlen2true%if we're also looking at 12.5 cm lengths in other versions
        newlinenum=linenum-itdif(1);
        if newlinenum<1
            newlinenum=linenum;
        end
        writematrix(EeX0end,filename,'Sheet',sheetnum,'Range',"AM"+num2str(newlinenum))        
        end
    else%otherwise all other irradiance vectors start at z
        writematrix(EeX0end,filename,'Sheet',sheetnum,'Range',"Z"+num2str(linenum))
    end
    %write everything stored for power accounting to out of the way place
    fibItname=fieldnames(FibIt);
    A=cell(1,length(fibItname));
    for w=1:(length(fibItname)-1)%
        A{1,w}=FibIt(jj).(fibItname{w})(iteration,h);
    end
    A{1,end}=cell2mat(FibIt(jj).Y(iteration,h));
        % struct2cell(FibIt(jj));%cellfun(@(f) FibIt(jj).(f)(iteration,h),FibIt(jj).(f), 'UniformOutput',false);%create an array of the lengths of each vector in each field
    
    writecell(A,filename,'Sheet',sheetnum,'Range',"CL"+num2str(linenum)) 

end
