classdef writeheadfile
    properties
        desc;%description, print to A1 description of whats happening on th
top;%'non iteratable parameters', print to B1

nonit;%values of each of the non iteratable parameters, print to C2

A6col;%'filename' 
B6col;%'iteration' title
iter;%iterable variable field titles, print to C4

reldata; %air sum side emission, transmission, water start side emission titles start this at O4
sim10;% titles to indicate start of simulated 10cm fiber trend, print to S3
xvec10;% x values of the 10cm simualtionm print to T4
sim50;% titles to indicate start of simulated 10cm fiber trend, print to AG3
xvec50;%x values of the 50cm fiber simulation print to AH4
xvecalt;%x values of other length fibers
simalt;%simulation of other fiber lengths
%remaining data fields will go after the 50cm variables:
energybal;%energy balance header title print to CF4
fieldFI;%tracked energy balance types print to CG4
    end
end