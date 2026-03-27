function  [iterParams,iterVars,paramNames,paramLengths,it_num,legendMain,Title_Main,yes_itname,it_name,itdif,xlen2true]= setupIterParams(iterParamsCustom,iterParamsStandard,SystemParam,comboit,readitlog,readitfile)


%get properties of the struct
iterParamsCustom.waterInterface=sort(iterParamsCustom.waterInterface,'ascend');
iterParamsCustom.xLen=sort(iterParamsCustom.xLen,'ascend');
iterParams = iterParamsCustom;%set = to the custom iteration struct. may change from here
paramNames = fieldnames(iterParams);%find the field names of the struct
paramLengths = cellfun(@(f) length(iterParams.(f)), paramNames);%create an array of the lengths of each vector in each field
paramNarray=string(paramNames);%convert the field names to a string array
iterVars = find(paramLengths > 1);
watlenvarind=length(paramNames);%should be the last variable
paramNames_=paramNames(1:end-1);%all of the param names except the water interface and xlen
iterVars_=iterVars(iterVars~=(watlenvarind));%all of the param iterations except water interface and xlen
itdif=zeros(2,1);
xlen2true=0;
watvar=1;
if length(iterParams.xLen)==2
    if iterParams.xLen(1)==12.5*10^4 && iterParams.xLen(2)==50.5*10^4
        iterVars_=iterVars_(iterVars_~=(watlenvarind-1));
        paramNames_=paramNames_(1:end-1);
        xlen2true=1;
    end
end
%if we're reading in the variables via a file
if readitlog==1%
    if ~isempty(readitfile)%if the file name is set
        filename=fullfile(pwd,'siminput',readitfile);%need the file to be created and stored in the siminput folder
        T=readtable(filename);%read the table from the file
        itvarnames=string(T.Properties.VariableNames);%iteration variable names
        skipf=zeros(1,length(itvarnames));%create a vector of the fields we should skip later
        ittot=skipf;%empty vector
        n=1;%index
        for k = itvarnames
            ittot(n)=length(T.(k));
            if iscell(T.(k))
                tempTk=zeros(ittot(n),1);
                for aa=1:ittot(n)

                   temptemp=str2num(T.(k){aa});
                   if length(temptemp)>1
                       tempv=sum(temptemp);
                   else
                       tempv=temptemp;
                   end
                   tempTk(aa)=tempv;%str2num(T.(k){aa});
                end
                iterParams.(k)=tempTk;
            else
            iterParams.(k)=T.(k);%store the table inpiut variable in the iterParams struct
            end
            skipf(n)=find(paramNarray==k);%store the string array index corresponding to thois field
            n=1+n;%increase index
        end
        if max(ittot)~=min(ittot)
            error('input file has different lengths of variables')
        end
        it_num=ittot(1);%all of the values should be equal, grab the first
        iterVars=skipf;
        %for the remaining parameters
        fieldnameL=1:length(paramNarray);
        fieldlogic=fieldnameL~=skipf(:);%logical index of params that havent been iterated over yet  #rows=length(skipf), #cols=length(fieldnameL)
        fields2write=paramNarray(sum(fieldlogic)==length(skipf));%only need to write the remaining variables (where the sum of each column to be equal to the #rows)
        for ll=1:length(fields2write)%for all the remaining fields
            %set the remaining iterable parameters to be a vector of the standard value, but with same length as the iteration
            m=fields2write(ll);
            iterParams.(m)=iterParamsStandard.(m).*ones(ittot(1),1);
        end
        it_name = ' ';
        yes_itname = 2;
        Title_Main = generateTitle(SystemParam);
         legendMain = '';
    else%the filename input is empty
        error('input simulation file name is empty')
        it_name = ' ';
        yes_itname = 0;
        legendMain = ' ';
        Title_Main = generateTitle(SystemParam);
        it_num=1;
    end
    xlen2true=0;
    if length(unique(iterParams.xLen))==2
    if iterParams.xLen(1)==12.5*10^4 && iterParams.xLen(2)==50.5*10^4
        xlen2true=1;
    end
    end
    if isscalar(unique(iterParams.waterInterface))
        watlenvarind=0;
    end
else
    % Check number of varying parameters
    if isempty(iterVars)
        it_name = ' ';
        yes_itname = 0;
        legendMain = ' ';
        Title_Main = generateTitle(SystemParam);
        it_num=1;
        %TODO: check the title function and set a description function
    elseif length(iterVars) > 1
        yes_itname = 2;
        if comboit==1%if we are iterating through all possible combinations

            % Multiple parameters vary — full grid
            [gridStruct, it_num] = generateParameterGrid(iterParamsCustom);
            legendMain = '';%generateLegends(gridStruct, paramNames(iterVars));
            
            it_name='';
            Title_Main = generateTitle(SystemParam);
            % Flatten grid to update iterParams values
            for i = 1:length(paramNames)
                name = paramNames{i};
                iterParams.(name) = gridStruct.(name);
            end
        else%we are iterating linearly
            %iterParams.(m)=iterParamsStandard.(m).*ones(ittot(1),1);
         
            it_num=sum(paramLengths(iterVars_))+1;
            startind=2;
            for i=1:length(paramNames_)
                assignvecbasic=iterParamsStandard.(paramNames_{i}).*ones(it_num,1);%set everything as a 1D vector length itnum of its standard value
                if any(i==iterVars_(:))
                    currentvar=iterParamsCustom.(paramNames_{i});
                    endind=startind+length(currentvar)-1;
                    assignvecbasic(startind:endind)=currentvar;%assign the custom params to the the vector
                    startind=endind+1;
                end
                if xlen2true
                    assignvecbasic=[assignvecbasic;assignvecbasic];%double it
                end
                if length(iterParams.waterInterface)>1
                    assignvecbasic=[assignvecbasic;assignvecbasic];%double it
                end
                iterParams.(paramNames_{i})=assignvecbasic;
            end
             if xlen2true
                it_num2=2*it_num;
                iterParams.xLen=(50.5*10^4).*ones(it_num2,1);
                iterParams.xLen(1:it_num,1)=(12.5*10^4).*ones(it_num,1);
                it_num=it_num2;
            end
            
            if length(iterParamsCustom.waterInterface)>1
                it_num2=2*it_num;
                iterParams.waterInterface=ones(it_num2,1);
                iterParams.waterInterface(1:it_num,1)=zeros(it_num,1);
                if length(iterParams.xLen)~=length(iterParams.waterInterface)
                    tempxlen=iterParams.xLen;
                    iterParams.xLen=[tempxlen;tempxlen];
                end
                it_num=it_num2;
            else
                iterParams.waterInterface=iterParamsStandard.waterInterface.*ones(it_num,1);
            end
          

            it_name = ' ';
            yes_itname = 2;
            legendMain = ' ';
            Title_Main = generateTitle(SystemParam);
        end
    else
        % Only one parameter is being iterated
        varName = paramNames{iterVars};
        values = iterParams.(varName);
        it_name = getDisplayName(varName);
        it_num = length(values);
        yes_itname = 1;
        Title_Main = generateTitle(SystemParam);
        legendMain = arrayfun(@(v) formatLegend(varName, v), values, 'UniformOutput', false);
        legendMain = string(legendMain);
    end
end
if xlen2true
    matIterVarCell=cell2mat(transpose(struct2cell(iterParams)));
    matIterVar=matIterVarCell(:,1:end-2);
    param1=matIterVar(1,:);
    param2=matIterVar(2,:);
    param1watlog=ismember(matIterVar,param1,'rows');
    param2watind=ismember(matIterVar,param2,'rows');
    indsP1=find(param1watlog==1);
    indsP2=find(param2watind==1);
    itdifxlen1=indsP1(2)-indsP1(1);
    itdifxlen2=indsP2(2)-indsP2(1);
    if itdifxlen1~=itdifxlen2
        error('weird xlen indexing issue 2')
    else
        itdif(1)=itdifxlen1;
    end
end


if length(iterParamsCustom.waterInterface)>1 && watlenvarind==1
    matIterVarCell=cell2mat(transpose(struct2cell(iterParams)));
    matIterVar=matIterVarCell(:,1:end-1);
    param1=matIterVar(1,:);
    param2=matIterVar(2,:);
    param1watlog=ismember(matIterVar,param1,'rows');
    param2watind=ismember(matIterVar,param2,'rows');
    indsP1=find(param1watlog==1);
    indsP2=find(param2watind==1);
    itdifairwat1=indsP1(2)-indsP1(1);
    itdifairwat2=indsP2(2)-indsP2(1);
    if length(indsP1)>2||length(indsP2)>2
        error('weird water indexing issue')
    if itdifairwat1~=itdifairwat2
        error('weird water indexing issue 2')
    else
        itdif(2)=itdifairwat1;
    end
end

end


