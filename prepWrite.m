function [header,topoffile] = prepWrite(filename,header,SystemParam,iterParam,FibIt,XVEC10,XVEC50,xvecall,desc,xlen2true)
%separate function to run out of the loop: write 2 file a table of non-iteratable variable settings

%if headercalc==0 %if we haven't calculated the header values yet then we need to calculate this
%start by defining all of the non-iteratable fields
fieldSys=string(fieldnames(SystemParam));
fieldIt=string(fieldnames(iterParam));
difLength=length(fieldSys)-length(fieldIt);
fieldRem=cell(1,difLength);%not sure if i need the +1 but it doesn't hurt to add
RemInd=zeros(1,difLength);
n=1;
for i=1:length(fieldSys)
    indfind=find(fieldIt==fieldSys(i));
    if ~any(indfind)%if there are no fieldIt matches with the current system param field
        fieldRem{n}=fieldSys(i);%store this as a non iteratable field parameter
        RemInd(n)=i;%store the index
        n=n+1;
    end
end
%get the data from those fields
RemVar=cell(size(RemInd));
for f=1:length(RemVar)
    % disp(f)
    % disp(fieldSys(RemInd(f)))
    % disp(SystemParam.(fieldSys(RemInd(f))))
     vartemp=SystemParam.(fieldSys(RemInd(f)));
    if ~isempty(vartemp)
        if ~iscell(vartemp)
            if imag(vartemp)==0
                RemVar{f}=vartemp;%SystemParam.(fieldSys(RemInd(f)));
            else
                RemVar{f}=string([num2str(real(vartemp)),'+ 1i.*',num2str(imag(vartemp),'%.2g')]);
            end
        else
            RemVar{f}='N/A';
        end
    else
        RemVar{f}='N/A';
    end    
end
sheetnum=1;%do this for the first sheet at least
%this is going to be specific to this particular simulation
header.desc=desc;%A1 description of whats happening on th
writematrix(header.desc,filename,'Sheet',sheetnum,'Range','A1');
header.top=['non iteratable parameters: ',string(fieldRem)];%print to B1
writematrix(header.top,filename,'Sheet',sheetnum,'Range','B1');
header.nonit=RemVar;%print to C2
writecell(header.nonit,filename,'Sheet',sheetnum,'Range','C2');
header.A6col='filename';
writematrix(header.A6col,filename,'Sheet',sheetnum,'Range','A6');
header.B6col='iteration';
writematrix(header.B6col,filename,'Sheet',sheetnum,'Range','B6');
header.iter=transpose(string(fieldnames(iterParam)));%print to C6
writematrix(header.iter,filename,'Sheet',sheetnum,'Range','C6');

header.reldata={'','Air Sum Side Emission (uW/cm2)','Transmitted (uW/cm2)','Water Sum Side Emission (uW/cm2)',''}; %start this at T6
writecell(header.reldata,filename,'Sheet',sheetnum,'Range','T6');
%% the XVEC print out up top will be dependent on if we're doing 12.5 and 50.5 cm or other lengths
if sum(xvecall==XVEC10)==length(xvecall) || sum(xvecall==XVEC50)==length(xvecall) || xlen2true
header.sim10={'10cm fiber sim Ee(x) uW/cm2';'X(cm)'};%print to Y5
header.xvec10=XVEC10;%print to Z6
writematrix(header.xvec10,filename,'Sheet',sheetnum,'Range','Z6');
header.sim50={'50cm fiber sim Ee(x) uW/cm2';'X(cm)'};%print to AL5
writecell(header.sim50,filename,'Sheet',sheetnum,'Range','AL5');
header.xvec50=XVEC50;%print to AM6
writematrix(header.xvec50,filename,'Sheet',sheetnum,'Range','AM6');
elseif isscalar(iterParms.xLen) 
    lenstr=sprintf('%d cm fiber sim Ee(x) uW/cm2',iterParams.xLen/1e4);
    header.simalt={lenstr;'X(cm)'};%print to Y5
    writecell(header.simalt,filename,'Sheet',sheetnum,'Range','Y5');
    writematrix(header.xvecalt,filename,'Sheet',sheetnum,'Range','Z6');

else
    lenstr=sprintf('variable length fiber (max %d cm)sim Ee(x) uW/cm2',max(iterParams.xLen/1e4));
    header.simalt={lenstr;'X(cm)'};%print to Y5
    writecell(header.simalt,filename,'Sheet',sheetnum,'Range','Y5');
    writematrix(header.xvecalt,filename,'Sheet',sheetnum,'Range','Z6');

end
%remaining data fields will go after the 50cm variables:
header.energybal='Energy balance data:';%print to CK6
writematrix(header.energybal,filename,'Sheet',sheetnum,'Range','CK6');
header.fieldFI=transpose(string(fieldnames(FibIt)));% print to CL6
writematrix(header.fieldFI,filename,'Sheet',sheetnum,'Range','CL6');


%create array/ table of iteratable variable settings
%write the header to the first sheet
topoffile=readcell(filename,'Sheet',1,'Range','A1:DC6');
end

%start the writing portion here
