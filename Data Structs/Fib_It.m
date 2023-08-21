classdef Fib_It %(Use FibIt(1).property for the avg, FibIt(2).property for the std dev)
    properties
        %power accounting
        Pow_enter ;%total power entering the optical fiber
        transmitted ;%total transmitted power
        pow_side ;% total side emitted power
        pow_side_use ;% total useable side emitted power
        SMA_pow_side ; %total side emitted power emitted w/in SMA
        pow_side_waterst;%total side emitted power within the water column
        %power accounting: losses
        absorbed ;%total UV absorbption loss
        backscat ;%total rayleigh loss
        SMAabs;% total power absorbed by SMA metal
        approxpow; %tracking the approximations i make. won't necessarily correspond to differences
        approxpow_dif;%sum of power differences due to difference in approximations
        approxpow_pos; %absolute value of all of the power difference approximations
        b2hpow;%total power transmitted back into the housing
        cutoffpow;%total power cut off to save model run time
        remaininglosses ;%total remaining losses of power unaccounted for
        %performance data
        UC;
        RatioIsIt;
        Y;
        
        
    end
end