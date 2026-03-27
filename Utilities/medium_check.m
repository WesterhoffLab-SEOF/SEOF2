function[n_medium,k_medium,transmission]=medium_check(SystemParam,x)
SMAFL=SystemParam.smaFlushLength;%=1*10^4;%the length of the SMA connector that is flush-ish to the fiber surface is 1 cm long

        %determine the refractive indices && attenuation coefficients
        if SystemParam.SMA==1 && x<=SMAFL && x>=0 %if it's within the flush length of the SMA connector, lose the light
            %select the refractive index and loss for the transmitting
            %material within the flush length of the SMA connector
            if SystemParam.isSmaSealed==1 && x <= SystemParam.smaFillLength %if we've sealed it with cytop, and are within the cytop fill region
                n_medium=SystemParam.n4;
                k_medium=SystemParam.kCytop;
            else %if it isn't sealed with cytop, or if we're in the part of the SMA "above" the cytop fill level
                n_medium=SystemParam.nMetal;
                k_medium=SystemParam.kAir;
            end
            transmission=0; %%either way, all the transmitted lightwill be lost to the
            %surface its touching, dont bother tracking "bounces"in the sma connector
        elseif x<0%this case shouldnt occur but if it does... we're in the internal housing
            n_medium=SystemParam.n2;
            k_medium=SystemParam.kAir;
            transmission=1;
        elseif SystemParam.isSmaSealed==1 && x <= SystemParam.smaFillLength %if we're within the sma connector fill portion and its sealed w cytop
            n_medium=SystemParam.n4;
            k_medium=SystemParam.kCytop;
            transmission=1;
        elseif SystemParam.waterInterface==1 && x < SystemParam.waterStart
            %if there is water, but we haven't reached it, then make sure it's
            %the air RI
            n_medium=SystemParam.n2;
            k_medium=SystemParam.kAir;
            transmission=1;
        else %we're in the selected main medium the fiber is submerged in
            n_medium=SystemParam.n5;
            k_medium=SystemParam.k;
            transmission=1;
        end
end