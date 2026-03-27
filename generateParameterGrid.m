function [gridStruct, it_num, iterVars] = generateParameterGrid(iterParams)
% Extract field names and values
paramNames = fieldnames(iterParams);
paramVals = (struct2cell(iterParams));

% Identify which parameters are iterable (i.e., length > 1)
iterFlags = cellfun(@(v) numel(v) > 1, paramVals);
iterVars = paramNames(iterFlags);
iterVals = paramVals(iterFlags);

% Handle three cases: 0, 1, or multiple iterable variables
nIterVars = numel(iterVals);

if nIterVars == 0
    % No iteration: just return single value struct
    gridStruct = iterParams;
    it_num = 1;
    return;

elseif nIterVars == 1
    % Single parameter to iterate
    v = iterVals{1};
    it_num = numel(v);
    gridStruct = repmat(struct(), 1, it_num);
    for i = 1:it_num
        for k = 1:numel(paramNames)
            pname = paramNames{k};
            if strcmp(pname, iterVars{1})
                gridStruct(i).(pname) = v(i);
            else
                gridStruct(i).(pname) = iterParams.(pname);
            end
        end
    end

else
    % Multiple parameters — use ndgrid
    meshh = combvec(paramVals{:});%ndgrid(cell2mat(paramVals))%iterVals{:})
    szmesh=size(meshh);
    it_num = szmesh(2);
    %gridStruct = repmat(struct(), 1, it_num)
    % for i = 1:it_num
    for k = 1:numel(paramNames)
        pname = paramNames{k};
        if ~any(imag(meshh(k,:)))
            paramit=real(meshh(k,:));
        else
            paramit=meshh(k,:);
        end
        gridStruct.(pname) = transpose(paramit);
    end
end
end

