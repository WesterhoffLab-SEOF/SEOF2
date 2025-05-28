function [FibIt, SumVal] = storeFiberResults(FibIt, SumVal, iteration, h, trial, numRandomIters)
    %storing average data 
    FibIt(1).Pow_enter(iteration,h)=mean(trial(1).Pow_enter(1,:));
    FibIt(1).transmitted(iteration,h)=mean(trial(1).transmitted(1,:));
    FibIt(1).pow_side(iteration,h)=mean(trial(1).pow_side(1,:));
    FibIt(1).pow_side_use(iteration,h)=mean(trial(1).pow_side_use(1,:));
    FibIt(1).SMA_pow_side(iteration,h)=mean(trial(1).SMA_pow_side(1,:));
    FibIt(1).pow_side_waterst(iteration,h)=mean(trial(1).pow_side_waterst(1,:));
    FibIt(1).absorbed(iteration,h)=mean(trial(1).absorbed(1,:));
    FibIt(1).backscat(iteration,h)=mean(trial(1).backscat(1,:));
    FibIt(1).SMAabs(iteration,h)=mean(trial(1).SMAabs(1,:));
    FibIt(1).approxpow(iteration,h)=mean(trial(1).approxpow(1,:));
    FibIt(1).approxpow_dif(iteration,h)=mean(trial(1).approxpow_dif(1,:));
    FibIt(1).approxpow_pos(iteration,h)=mean(trial(1).approxpow_pos(1,:));
    FibIt(1).b2hpow(iteration,h)=mean(trial(1).b2hpow(1,:));
    FibIt(1).cutoffpow(iteration,h)=mean(trial(1).cutoffpow(1,:));
    FibIt(1).remaininglosses(iteration,h)=mean(trial(1).remaininglosses);
    
    %std deviations
    FibIt(2).Pow_enter(iteration,h)=std(trial(1).Pow_enter(1,:));
    FibIt(2).transmitted(iteration,h)=std(trial(1).transmitted(1,:));
    FibIt(2).pow_side(iteration,h)=std(trial(1).pow_side(1,:));
    FibIt(2).pow_side_use(iteration,h)=std(trial(1).pow_side_use(1,:));
    FibIt(2).SMA_pow_side(iteration,h)=std(trial(1).SMA_pow_side(1,:));
    FibIt(2).pow_side_waterst(iteration,h)=std(trial(1).pow_side_waterst(1,:));
    FibIt(2).absorbed(iteration,h)=std(trial(1).absorbed(1,:));
    FibIt(2).backscat(iteration,h)=std(trial(1).backscat(1,:));
    FibIt(2).SMAabs(iteration,h)=std(trial(1).SMAabs(1,:));
    FibIt(2).approxpow(iteration,h)=std(trial(1).approxpow(1,:));
    FibIt(2).approxpow_dif(iteration,h)=std(trial(1).approxpow_dif(1,:));
    FibIt(2).approxpow_pos(iteration,h)=std(trial(1).approxpow_pos(1,:));
    FibIt(2).b2hpow(iteration,h)=std(trial(1).b2hpow(1,:));
    FibIt(2).cutoffpow(iteration,h)=std(trial(1).cutoffpow(1,:));
    FibIt(2).remaininglosses(iteration,h)=std(trial(1).remaininglosses);

    %average metrics
    if numRandomIters==1
         FibIt(1).Y(iteration,h)=trial(1).Y;
    else
        FibIt(1).Y(iteration,h)={mean(cell2mat(trial(1).Y),2)};
        FibIt(1).Y(iteration,h)={mean(cell2mat(trial(1).Y),2)};
    end


    FibIt(1).UC(iteration,h)=mean(trial(1).UC(1,:));
    FibIt(1).RatioIsIt(iteration,h)=mean(trial(1).RatioIsIt(1,:));
    FibIt(2).UC(iteration,h)=std(trial(1).UC(1,:));
    FibIt(2).RatioIsIt(iteration,h)=std(trial(1).RatioIsIt(1,:));

    % % Storing the larger data for maybe checking the trends out with
    % SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
    % SumVal.Pow_enter(iteration,h)={ray.Pow_enter};
    % SumVal.transmitted(iteration,h)={ray.transmitted};
    % SumVal.pow_side(iteration,h)={ray.pow_side};
    % SumVal.pow_side_use(iteration,h)={ray.pow_side_use};
    % SumVal.SMA_pow_side(iteration,h)={ray.SMA_pow_side};
    % SumVal.pow_side_waterst(iteration,h)={ray.pow_side_waterst};
    % SumVal.absorbed(iteration,h)={ray.absorbed};
    % SumVal.backscat(iteration,h)={ray.backscat};
    % SumVal.SMAabs(iteration,h)={ray.SMAabs};
    % SumVal.approxpow_dif(iteration,h)={ray.approxpow_dif};%sum of power differences due to difference in approximations
    % SumVal.approxpow_pos(iteration,h)={ray.approxpow_pos}; %absolute value of all of the power difference approximationsSumVal.b2hpow(iteration,h)={ray.;
    % SumVal.cutoffpow(iteration,h)={ray.cutoffpow};
    % SumVal.remaininglosses(iteration,h)={trial(1).remaininglosses};
    % SumVal.UC(iteration,h)={ray.UC};
    % SumVal.RatioIsIt(iteration,h)={ray.RatioIsIt};
    % SumVal.Y(iteration,h)={ray.Y};
         

end