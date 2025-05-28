function [ray, trial] = initializeTrial(a, b, numRandomIters)
    trial(1).Pow_enter=zeros(1,numRandomIters);
    trial(1).transmitted=zeros(1,numRandomIters);
    trial(1).pow_side=zeros(1,numRandomIters);
    trial(1).pow_side_use=zeros(1,numRandomIters);
    trial(1).SMA_pow_side=zeros(1,numRandomIters);
    trial(1).pow_side_waterst=zeros(1,numRandomIters);
    trial(1).absorbed=zeros(1,numRandomIters);
    trial(1).backscat=zeros(1,numRandomIters);
    trial(1).SMAabs=zeros(1,numRandomIters);
    trial(1).approxpow=zeros(1,numRandomIters);
    trial(1).approxpow_dif=zeros(1,numRandomIters);
    trial(1).approxpow_pos=zeros(1,numRandomIters);
    trial(1).b2hpow=zeros(1,numRandomIters);
    trial(1).cutoffpow=zeros(1,numRandomIters);
    trial(1).remaininglosses=zeros(1,numRandomIters);
    trial(1).UC=zeros(1,numRandomIters);trial(2).UC=zeros(1,numRandomIters);
    trial(1).RatioIsIt=zeros(1,numRandomIters);trial(2).RatioIsIt=zeros(1,numRandomIters);
    trial(1).Y=cell(numRandomIters,1);trial(2).Y=cell(numRandomIters,1);
    % trial(1).meas_inten=cell(1,numRandomIters); %length of all the rays used
    % trial(1).meas_points=cell(1,numRandomIters);%length of all the rays used
    trial(1).meas_inten=cell((a*b),numRandomIters); %length of all the rays used
    trial(1).meas_points=cell((a*b),numRandomIters);%length of all the rays used
    %empty storage vector setup for each individual ray 
    %possibly will need to include dimension for the scatter cone later
    %(numRandomIters,a,b,gg_max)
    ray.Pow_enter=zeros(numRandomIters,a,b);%=I_ent0
    ray.transmitted=zeros(numRandomIters,a,b);
    ray.pow_side=zeros(numRandomIters,a,b);
    ray.pow_side_use=zeros(numRandomIters,a,b);
    ray.SMA_pow_side=zeros(numRandomIters,a,b);
    ray.pow_side_waterst=zeros(numRandomIters,a,b);
    ray.absorbed=zeros(numRandomIters,a,b);
    ray.backscat=zeros(numRandomIters,a,b);
    ray.SMAabs=zeros(numRandomIters,a,b);
    ray.approxpow=zeros(numRandomIters,a,b);
    ray.approxpow_dif=zeros(numRandomIters,a,b);%sum of power differences due to difference in approximations
    ray.approxpow_pos=zeros(numRandomIters,a,b); %absolute value of all of the power difference approximations
    ray.b2hpow=zeros(numRandomIters,a,b);
    ray.cutoffpow=zeros(numRandomIters,a,b);
    ray.remaininglosses=zeros(numRandomIters,a,b);%should theoretically be equal to the approxpow_dif?;
    ray.UC=zeros(numRandomIters,a,b);
    ray.RatioIsIt=zeros(numRandomIters,a,b);
    ray.Y=cell(numRandomIters,a,b);
    ray.meas_inten=cell(numRandomIters,a,b);
    ray.meas_points=cell(numRandomIters,a,b);
end