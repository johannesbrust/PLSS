% All methods for overdetermined

function OUTPUTS= all_overdet_methods(OUTPUTS, options, Prob,num_methods,max_iterations)

if(num_methods <1) return; end
%% Gaussian Least squares
options.max_iterations =max_iterations(1);
[x, output_gaussls] = solve_system(Prob,@iter_Gaussls, @boot_Gaussls,options );
OUTPUTS = [ OUTPUTS ; output_gaussls];
if(num_methods <2) return; end
%% Gaussian Iterative Projection IGP
options.max_iterations =max_iterations(2);
[x, output_gaussigp] = solve_system(Prob,@iter_Gaussigp, @boot_Gaussigp,options );
OUTPUTS = [ OUTPUTS ; output_gaussigp];
if(num_methods <3) return; end
%% LS Coordinate Descent
options.max_iterations =max_iterations(3);
[x, output_lscd] = solve_system(Prob,@iter_LSCD, @boot_LSCD,options );
OUTPUTS = [ OUTPUTS ; output_lscd];
if(num_methods <4) return; end
%% Randomized Karczmarz
options.max_iterations =max_iterations(4);
 [x, output_kz] = solve_system(Prob,@iter_Karczmarz, @boot_Karczmarz,options );
 OUTPUTS = [ OUTPUTS ; output_kz];
end