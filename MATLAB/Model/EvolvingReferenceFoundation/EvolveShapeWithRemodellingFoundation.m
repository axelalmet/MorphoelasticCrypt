function EvolveShapeWithRemodellingFoundation
% Set the parameters
kf = 0.16; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 2*sqrt(3)*L0/h; %Dimensionless length
y0 = 6*L;
K = kf*h/(12*w); % Dimensionless foundation stiffness
Es = 1; % Stretching stiffness
b1 = 0.5; % Bending stiffness
dt = 1e-3; % Time step

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsk0p02L29_sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

solFromData.y(3,:) = y0 + solFromData.y(3,:);


sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% % Define the Wnt function
% W = @(S, width) exp(-((S - 0.5*L)/width).^2);
W = @(S, mu1, mu2, width) 0.5.*(tanh((S - mu1)./width) - tanh((S - mu2)./width));

% Set growth parameters
mu1 = 0.5*L - 2*sqrt(3)*w/h;
mu2 = 0.5*L + 2*sqrt(3)*w/h;
width = 0.1;
sigma = mu2 - mu1;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

parameters.mu1 = mu1;
parameters.mu2 = mu2;
parameters.width = 0.1;
parameters.currentArcLength = solFromData.y(1,:);
parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.eta = eta; % Rate of chemical change
parameters.Es = Es; % Stretch stiffness
parameters.b1 = b1;
parameters.Eb = 1; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu = etaF*eta; % Foundation relaxation timescale
parameters.etaK = kf/(eta*etaF); % Curvature relaxation timescale
parameters.dt = dt; % Time step

%% Solve the initial bvp to obtain a structure for the first solution.
SOld = solFromData.y(1,:);

% Initialise non-uniform gamma
gammaOld = 1;
% firstGamma = gammaOld.*(1 + dt*W(solFromData.y(1,:), sigma));
firstGamma = gammaOld.*(1 + dt.*W(solFromData.y(1,:), mu1, mu2, width));
parameters.gamma = firstGamma;

parameters.Eb = 1 - b1.*exp(-((solFromData.y(1,:) - 0.5*L)./sigma).^2);

% Initialise foundation shape
parameters.Px = SOld;
parameters.Py = y0.*ones(1, length(solFromData.x));

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters);

tic
% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);

toc
% plot(initSol.x, initSol.y(3,:))

initSol.x = solFromData.x;
initSol.y = deval(numSol, solFromData.x); 
                            
%%  
% Update the initial data
solMesh = solFromData.x;

gammaOld = interp1(solFromData.x, firstGamma, initSol.x);
parameters.gamma = gammaOld;
parameters.currentArcLength = interp1(solFromData.x, cumtrapz(solFromData.y(1,:), parameters.gamma), initSol.x);
parameters.mu1 = 0.5*(parameters.currentArcLength(end)) - 2*sqrt(3)*w/h;
parameters.mu2 = 0.5*(parameters.currentArcLength(end)) + 2*sqrt(3)*w/h;
parameters.sigma = (parameters.mu2) - (parameters.mu1);
parameters.Eb =  1 - b1.*exp(-((parameters.currentArcLength - 0.5*(parameters.mu1 + parameters.mu2))./parameters.sigma).^2);

%%

solOld = initSol;

initX = solOld.y(1,:);
initY = solOld.y(2,:);

parameters.Px = (1 - dt*(parameters.nu)).*SOld + dt*(parameters.nu).*initX;
parameters.Py = (1 - dt*(parameters.nu)).*(parameters.y0) + dt*(parameters.nu).*initY;
    
% Set the times we want to solve the problem for
dt = 2.5*1e-2;
parameters.dt = dt; 

TMax = 5.0;
times = [0, 1e-3:dt:TMax];
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);
gammaSols = cell(numSols, 1);
foundationSols = cell(numSols, 1);
uHatSols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;  

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
uHatSols{1} = [flatSol.x; 0.*flatSol.x];
foundationSols{1} = [initSol.y(1,:); y0.*ones(1, length(initSol.x))];


% First non-trivial solution
Sols{2} = [initSol.x; initSol.y];
gammaSols{2} = [L.*initSol.x; gammaOld];
uHatSols{2} = [L.*initSol.x, parameters.uHat];
foundationSols{2} = [parameters.Px; parameters.Py];

tic

% Update the solutions in time
for i = 3:numSols
        
    % Update the solution
    [solMeshNew, solNew, gammaNew, EbNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolution(solMesh, solOld, W, parameters, solOptions);
    
    % Update the solutions, gamma, and the spring stresses
    gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    parameters.gamma = gammaOld;
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    parameters.Eb = interp1(solOld.x, EbNew, solNew.x);
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma);
    
    parameters.mu1 = 0.5*(parameters.currentArcLength(end)) - 2*sqrt(3)*w/h;
    parameters.mu2 = 0.5*(parameters.currentArcLength(end)) + 2*sqrt(3)*w/h; 
    parameters.sigma = parameters.mu2 - parameters.mu1;
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (trapz(solNew.y(1,:), gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
        
        Sols = Sols(1:(i - 1));
        gammaSols = gammaSols(1:(i - 1));
        times = times(1:(i - 1));
        foundationSols = foundationSols(1:(i - 1));
        uHatSols = uHatSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;

    Sols{i} = [solOld.x; solOld.y];
    gammaSols{i} = [L.*solNew.x; gammaOld];
    foundationSols{i} = [parameters.Px; parameters.Py];
    uHatSols{i} = [L.*solNew.x; parameters.uHat];
    
    solMesh = solMeshNew;
                        
    end

toc

% Save the solutions
outputDirectory = '../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_current_0p5_sigmaE_2w_nu_1_k_0p02_L0_0p125_current_step_sigma_0p1_init_2w_etaK_0p16';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
