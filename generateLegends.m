function legends = generateLegends(gridStruct, fieldnames)
    n = numel(gridStruct);  
    legends = strings(1, n);
    for i = 1:n
        parts = cellfun(@(f) formatLegend(f, gridStruct(i).(f)), fieldnames, 'UniformOutput', false);
        legends(i) = strjoin(parts, ', ');
    end
end
