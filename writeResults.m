function res = writeResults(SystemParam, description, filename, sheet_num, iteration, h, FibIt, c, legend_Main, aa_lim, LED_d, XVEC)

if (SystemParam.waterInterface == 0)
    ext_media="External Medium: Air"; %string for excel document description
else
    ext_media="External Medium: Water"; %string for excel document description
end

%documentation of variables used in the 1st row
if c>1 %if multiple fibers
    description=description + " fiber no. " +h;
end
writematrix(description, filename,'Sheet', sheet_num + (h-1),'Range', "A1");
writematrix(" iteration; " + iteration +", "+ legend_Main(iteration), filename,'Sheet', sheet_num + (h-1),'Range', "A"+ num2str((6*iteration)-2));
writematrix("Separation Distance: " + SystemParam.led_distance + " (um)", filename,'Sheet', sheet_num + (h-1),'Range', "B1");
writematrix("Wavelength; " + SystemParam.uv_wavelength + "(nm)", filename,'Sheet', sheet_num + (h-1),'Range', "C1");
writematrix("Angle resolution: " + SystemParam.angnum, filename,'Sheet', sheet_num + (h-1),'Range', "E1");
writematrix("initial intensity Entering the Fiber: " + SystemParam.I_init + " (mW)", filename,'Sheet', sheet_num + (h-1),'Range', "F1");
%         writematrix("Polymer coat: " + SystemParam.poly_coat, filename,'Sheet', sheet_num + (h-1),'Range', "G1");
%         writematrix("Polymer thickness (if true): " + polymer_thickness + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "H"+ num2str(iteration));
%         writematrix("Nanoparticle Density: " + np_density, filename,'Sheet', sheet_num + (h-1),'Range', "I1");
%         writematrix("Nanoparticle Layer Thickness: " + np_thickness + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "J1");
%         writematrix("Nanoparticle Diameter: " + SystemParam.particle_dmtr + " (nm)", filename,'Sheet', sheet_num + (h-1),'Range', "K1");
%writematrix("# of Mie vectors considered: " + num_vectors, filename,'Sheet', sheet_num + (h-1),'Range', "L1");   
%writematrix("Measurement distance: " + SystemParam.meas_distance + ( "um"), filename,'Sheet', sheet_num + (h-1),'Range', "M1");
writematrix(ext_media, filename,'Sheet', sheet_num + (h-1),'Range', "N1");
writematrix("# of for loops used to smooth random noise: " + aa_lim, filename,'Sheet', sheet_num + (h-1),'Range', "O1");
writematrix("Length of fiber division in plot: " + (SystemParam.division/10e3) + " (cm)", filename,'Sheet', sheet_num + (h-1),'Range', "P1");
writematrix("Photons Considered: " +SystemParam.raynum + ", Angles Considered: "+  SystemParam.angnum, filename,'Sheet', sheet_num + (h-1),'Range', "Q1");
writematrix("Number Fibers: " +SystemParam.num_fib, filename,'Sheet', sheet_num + (h-1),'Range', "R1");

housepow=1;
%iteration values

writematrix("Fiber length: " + (SystemParam.xlen/10000) + " (cm)", filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)-3));
writematrix("Radius of fiber: " + SystemParam.r_fiber + " (um)", filename,'Sheet', sheet_num + (h-1),'Range', "B" + num2str((6*iteration)-3));
writematrix("Distance of LED from Fiber: " +LED_d + " (um)", filename,'Sheet', sheet_num + (h-1),'Range', "C"+ num2str((6*iteration)-3));
writematrix("Complex Refractive Index of Fiber:" + SystemParam.n1, filename,'Sheet', sheet_num + (h-1),'Range', "D"+ num2str((6*iteration)-3));
writematrix("Complex Refractive Index of Metal:" + SystemParam.n_metal, filename,'Sheet', sheet_num + (h-1),'Range', "D"+ num2str((6*iteration)-2));

%documentation of the data generated in the iteration
writematrix("Total light entering the fiber: " + FibIt(1).Pow_enter(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "E" + num2str((6*iteration)-3));
writematrix("Total transmitted light; " + FibIt(1).transmitted(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "F" + num2str((6*iteration)-3));
writematrix("Total side emitted light; " + FibIt(1).pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "G" + num2str((6*iteration)-3));
writematrix("Total useable side emitted light; " + FibIt(1).pow_side_use(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "H" + num2str((6*iteration)-3));
writematrix("Total side emitted light w/in the SMA: " + FibIt(1).SMA_pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "I" + num2str((6*iteration)-3));
writematrix("Total adsorbed light: " + FibIt(1).absorbed(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "J" + num2str((6*iteration)-3));
writematrix("Total back scattered light: " + FibIt(1).backscat(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "K" + num2str((6*iteration)-3));
writematrix("Total SMA adsorbed light: " + FibIt(1).SMAabs(iteration,h)+ " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "L" + num2str((6*iteration)-3));
writematrix("Total approximation loss light: " + FibIt(1).approxpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "M" + num2str((6*iteration)-3));
writematrix("Total back 2housing loss light: " + FibIt(1).b2hpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "N" + num2str((6*iteration)-3));
writematrix("Total cutoff loss light: " + FibIt(1).cutoffpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "O" + num2str((6*iteration)-3));
writematrix("Total remaining light loss: " + FibIt(1).remaininglosses(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "P" + num2str((6*iteration)-3));
writematrix("Avg Uniformity Cofficient: " + FibIt(1).UC(iteration,h), filename,'Sheet', sheet_num + (h-1),'Range', "Q" + num2str((6*iteration)-3));
writematrix("Avg Is/It Ratio: " + FibIt(1).RatioIsIt(iteration,h) , filename,'Sheet', sheet_num + (h-1),'Range', "R" + num2str((6*iteration)-3));
writematrix("Total water st side emitted light; " + FibIt(1).pow_side_waterst(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "S" + num2str((6*iteration)-3));

%documentation of the std deviation of the generated data
if aa_lim>1
    writematrix(FibIt(2).Pow_enter(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "E" + num2str((6*iteration)-2));
    writematrix(FibIt(2).transmitted(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "F" + num2str((6*iteration)-2));
    writematrix(FibIt(2).pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "G" + num2str((6*iteration)-2));
    writematrix(FibIt(2).pow_side_use(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "H" + num2str((6*iteration)-2));
    writematrix(FibIt(2).SMA_pow_side(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "I" + num2str((6*iteration)-2));
    writematrix(FibIt(2).absorbed(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "J" + num2str((6*iteration)-2));
    writematrix(FibIt(2).backscat(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "K" + num2str((6*iteration)-2));
    writematrix(FibIt(2).SMAabs(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "L" + num2str((6*iteration)-2));
    writematrix(FibIt(2).approxpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "M" + num2str((6*iteration)-2));
    writematrix(FibIt(2).b2hpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "N" + num2str((6*iteration)-2));
    writematrix(FibIt(2).cutoffpow(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "O" + num2str((6*iteration)-2));
    writematrix(FibIt(2).remaininglosses(iteration,h) + " (uW)", filename,'Sheet', sheet_num + (h-1),'Range', "P" + num2str((6*iteration)-2));
    writematrix(FibIt(2).UC(iteration,h), filename,'Sheet', sheet_num + (h-1),'Range', "Q" + num2str((6*iteration)-2));
    writematrix(FibIt(2).RatioIsIt(iteration,h) , filename,'Sheet', sheet_num + (h-1),'Range', "R" + num2str((6*iteration)-2));
    writematrix(FibIt(2).pow_side_waterst(iteration,h), filename,'Sheet', sheet_num + (h-1),'Range', "S" + num2str((6*iteration)-2));
end
%side emission over a the fiber distance and iteration
writematrix(XVEC,filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)-1));
writematrix(cell2mat(FibIt(1).Y(iteration,h)),filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)));
writematrix(cell2mat(FibIt(2).Y(iteration,h)),filename,'Sheet', sheet_num + (h-1),'Range', "A" + num2str((6*iteration)+1));