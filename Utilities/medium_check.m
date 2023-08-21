function[n_medium,k_medium]=medium_check(SystemParam,x)
        %determine the refractive indices && attenuation coefficients
        if SystemParam.waterInterface==1 && x < SystemParam.waterstart
            %if there is water, but we haven't reached it, then make sure it's
            %the air RI
            n_medium=SystemParam.n2;
            k_medium=SystemParam.kair;
        else%otherwise the RI/attenuation coeff is that of the medium
            n_medium=SystemParam.n5;
            k_medium=SystemParam.k;%
        end
end