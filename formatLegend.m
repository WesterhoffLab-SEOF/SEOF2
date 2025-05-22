function str = formatLegend(param, val)
    if islogical(val)
        str = sprintf('%s: %s', param, mat2str(val));
    elseif isreal(val)
        str = sprintf('%s: %.3g', param, val);
    else
        str = sprintf('%s: %.3g + %.3gi', param, real(val), imag(val));
    end
end