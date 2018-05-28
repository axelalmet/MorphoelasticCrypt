function EvolveLinearKelvinVoigtStressesWithoutGrowth
% Load the solutions previously computed
outputValues = 'Eb_0p75_sigmaE_2w_nu_4p8_k_0p02_L0_0p125_sigma_2w_area_1_mu_0_inext_currentforce';
outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/KelvinVoigt/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'kvstresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses
initSol.x = Sols{end - 2}(1,:);
initSol.y = Sols{end - 2}(2:end,:);
initGamma = gammaSols{end - 2};
initP = stressSols{end - 2};

parameters.gamma = initGamma(2,:);
parameters.P = initP(2,:);

solMesh = initSol.x;

% Set the solver options for bvp4c
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

%% Update the stresses until self-intersection or they tend to a steady state

TMax = 4.5; 
dt = parameters.dt;
L = parameters.L;
newTimes = [0, 1e-3:dt:TMax];
numSols = length(newTimes);


SolsWithoutGrowth = cell(numSols, 1);
stressSolsWithoutGrowth = cell(numSols, 1);

solOld = initSol;

SolsWithoutGrowth{1} = Sols{end};
stressSolsWithoutGrowth{1} = stressSols{end};

tic

for i = 2:numSols
        
    % Update the solution
    [solNew, PNew] = UpdateLinearKelvinVoigtSolutionWithoutGrowth(solMesh, solOld, parameters, solOptions);     
        
    % Stop the solution the curve self-intersects or the stress tends to a
    % steady state
    if ( (~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))||(norm(PNew - parameters.P) < 1e-6) )
        
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        stressSolsWithoutGrowth = stressSolsWithoutGrowth(1:(i - 1));
        
        break
    end
    
    % Update the solutions and the spring stresses
    parameters.P = interp1(solOld.x, PNew, solNew.x);
    solOld = solNew;

    SolsWithoutGrowth{i} = [solOld.x; solOld.y];
    stressSolsWithoutGrowth{i} = [L.*solNew.x; parameters.P];
                        
end

toc

outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/KelvinVoigt/'; 
outputValues = 'Eb_0p75_sigmaE_2w_nu_4p8_k_0p02_L0_0p125_sigma_2w_area_1_mu_0_inext_currentforce_nogrowth';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
save([outputDirectory, 'kvstresses_', outputValues,'.mat'], 'stressSolsWithoutGrowth') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times