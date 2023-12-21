function[n_medium,k_medium,transmission]=medium_check(SystemParam,x)
SMAFL=SystemParam.SMA_flushlength;%=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long

        %determine the refractive indices && attenuation coefficients
        if SystemParam.SMA==1 && x<=SMAFL && x>=0 %if it's within the flush length of the SMA connector, lose the light
        
            n_medium=SystemParam.n_metal;
            k_medium=SystemParam.kair;
            transmission=0; %multiply the transmission by this
        elseif x<0%this case shouldnt occur but if it does... it's the internal housing
            n_medium=SystemParam.n2;
            k_medium=SystemParam.kair;
            transmission=1;
        elseif SystemParam.sealed_SMA==1 && x <= SystemParam.SMA_fill %if we're within the sma connector fill portion and its sealed w cytop
            n_medium=SystemParam.n4;
            k_medium=SystemParam.k_cytop;
            transmission=1;
        elseif SystemParam.waterInterface==1 && x < SystemParam.waterstart
            %if there is water, but we haven't reached it, then make sure it's
            %the air RI
            n_medium=SystemParam.n2;
            k_medium=SystemParam.kair;
            transmission=1;
        else
%         elseif x <= SystemParam.xlen%otherwise the RI/attenuation coeff is that of the medium
            n_medium=SystemParam.n5;
            k_medium=SystemParam.k;%
            transmission=1;
%         else
%             error('x out of bounds')
        end
end
