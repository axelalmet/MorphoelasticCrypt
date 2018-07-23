function EvolveShapeWithLinearMaxwellFoundation
% Set the parameters
kf = 0.16; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 2*sqrt(3)*L0/h; %Dimensionless length
y0  = 6*L;
K = kf*h/(12*w); % Dimensionless foundation stiffness
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
b1 = 0.5; % Bending stiffness
dt = 1e-3; % Time step

% Get the initial solution from AUTO 
solData = load('../../../../Data/planarmorphorodsk0p02L29_sol_1'); %
% Cartesian basis

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

% Translate solution up
solFromData.y(3,:) = y0 + solFromData.y(3,:);

% 
sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% % Define the Wnt function
W = @(S, width) exp(-((S - 0.5*L)./width).^2);

eta = 1/24; % Define eta such that \eta^{-1} = 24 hours
mu = 0; 
nu = eta^(-1);

parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.mu = mu; % Rate of mechanical inhibition
parameters.eta = eta; % Rate of chemical change
parameters.n3s = n3s; % Target axial stress
parameters.Es = Es; % Stretch stiffness
parameters.b1 = b1;
parameters.Eb = 1 - b1.*W(solFromData.y(1,:), sigma); % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu = kf/(eta*nu); % Foundation relaxation timescale
parameters.etaK = parameters.nu; % Curvature relaxation timescale
parameters.dt = dt; % Time step

%% Solve the initial bvp to obtain a structure for the first solution.

% Cartesian basis solution extraction
yOld = solFromData.y(3,:);
FOld = solFromData.y(4,:);
GOld = solFromData.y(5,:);
thetaOld = solFromData.y(6,:);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

gammaOld = 1;
% firstGamma = gammaOld.*(1 + dt*(W(solFromData.y(1,:), sigma./gammaOld) + mu.*(n3Old - n3s)));
firstGamma = gammaOld.*(1 + dt*(W(solFromData.y(1,:), sigma)));
parameters.gamma = firstGamma;

parameters.Eb = 1 - b1.*W(solFromData.y(1,:), sigma./firstGamma);
% parameters.Eb = 1 - b1.*W(solFromData.y(1,:), sigma);

parameters.P = yOld - (parameters.y0);
parameters.uHat = zeros(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) SimplifiedLinearMaxwellFoundationOdes(x, M, solFromData, parameters);

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
gammaOld = interp1(solFromData.x, firstGamma, initSol.x);
parameters.gamma = gammaOld;
solOld = initSol;
    
initY = initSol.y(3,:);

parameters.P = initY - (parameters.y0);
parameters.uHat = zeros(1, length(solFromData.x));

solMesh = solFromData.x;
    
% Set the times we want to solve the problem for
dt = 2.5*1e-2;
parameters.dt = dt; 

TMax = 4.5;
times = [0, 1e-3:dt:TMax];
numSols = length(times);    

% Initialise the solutions
Sols = cell(numSols, 1);    
gammaSols = cell(numSols, 1);
stressSols = cell(numSols, 1);
uHatSols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;      

flatSol.y(3,:) = 0.*flatSol.y(3,:) + y0;
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
stressSols{1} = [flatSol.x; flatSol.y];
uHatSols{1} = [flatSol.x; flatSol.y];


% First non-trivial solution
Sols{2} = [initSol.x; initSol.y];
gammaSols{2} = [L.*initSol.x; gammaOld];
stressSols{2} = [L.*initSol.x; parameters.P];
uHatSols{2} = [L.*initSol.x, parameters.uHat];

tic

% Update the solutions in time
for i = 3:numSols
    
    % Update the solution
    [solMeshNew, solNew, gammaNew, EbNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellSolution(solMesh, solOld, W, parameters, solOptions);
    
    % Update the solutions, gamma, and the spring stresses
    gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    parameters.gamma = gammaOld;
    parameters.P = interp1(solOld.x, PNew, solNew.x);
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    parameters.Eb = interp1(solOld.x, EbNew, solNew.x); 
        
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (trapz(solNew.x, gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
        
        Sols = Sols(1:(i - 1));
        gammaSols = gammaSols(1:(i - 1));
        times = times(1:(i - 1));
        stressSols = stressSols(1:(i - 1));
        uHatSols = uHatSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    
    Sols{i} = [solOld.x; solOld.y];
    gammaSols{i} = [L.*solNew.x; gammaOld];
    stressSols{i} = [L.*solNew.x; parameters.P];
    uHatSols{i} = [L.*solNew.x; parameters.uHat];
    
    solMesh = solMeshNew;
    
end

toc

% Save the solutions

outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/Maxwell/';
outputValues = 'Eb_current_0p5_sigmaE_2w_nu_simple_init_0p16_k_0p02_L0_0p125_sigma_init_2w_etaK_0p16';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'maxwellstresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

