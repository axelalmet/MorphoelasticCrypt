function EvolveShapeWithHybridFoundation
% Set the parameters
kf = 0.16; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 2*sqrt(3)*L0/h; %Dimensionless length
K = kf*h/(12*w); % Dimensionless foundation stiffness
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
Eb = 1.0; % Bending stiffness
dt = 1e-4; % Time step

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsdirectorbasisk0p02L29_sol_3');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

% Director basis solution extraction
SOld = solFromData.y(1,:);
r1Old = solFromData.y(2,:);
r3Old = solFromData.y(3,:);
n3Old = solFromData.y(5,:);
thetaOld = solFromData.y(6,:);

xOld = -r1Old.*sin(thetaOld) + r3Old.*cos(thetaOld);
yOld = r1Old.*cos(thetaOld) + r3Old.*sin(thetaOld);

% 
sigma = 0.1*L; % "Width" of Wnt gradient
% % Define the Wnt function
W = @(S, width) exp(-(L*(S - 0.5)/width).^2);

eta = 1.0/trapz(solFromData.x, W(solFromData.x, sigma)); % Define eta such that the area is unit one
% eta = 1;
mu = 0; 
nu1 = eta^(-1);
nu3 = eta^(-1);

parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.sigma = sigma; % Width of wnt gradient
parameters.mu = mu; % Rate of mechanical inhibition
parameters.eta = eta; % Rate of chemical change
parameters.n3s = n3s; % Target axial stress
parameters.Es = Es; % Stretch stiffness
parameters.Eb = Eb; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu1 = K/(eta*nu1); % Foundation relaxation timescale
parameters.nu3 = K/(eta*nu3); % Foundation relaxation timescale
parameters.dt = dt; % Time step

%% Solve the initial bvp to obtain a structure for the first solution.

gammaOld = 1;
firstGamma = gammaOld.*(1 + dt*(W(solFromData.x, sigma) + mu.*(n3Old - n3s)));
parameters.gamma = firstGamma;

% parameters.K = K.*firstGamma;

% parameters.P1 = dt.*(parameters.nu1).*r1Old;
% parameters.P3 = r3Old - SOld.*cos(thetaOld);

% Define the ODEs and BCs
DerivFun = @(x, M) FlatFoundationNonlinearGeomOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) NonlinearGeomBCs(Ml, Mr, parameters);

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
parameters.K = K.*gammaOld;
solOld = initSol;

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

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;  

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];

% First non-trivial solution
Sols{2} = [initSol.x; initSol.y];
gammaSols{2} = [L.*initSol.x; gammaOld];

tic

% Update the solutions in time
for i = 3:numSols
        
    % Update the solution
    [solNew, gammaNew, PNew] = UpdateLinearViscoelasticSolution(solMesh, solOld, W, parameters, solOptions);     
    
    % Update the solutions, gamma, and the spring stresses
    gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    parameters.gamma = gammaOld;
    parameters.PX = interp1(solOld.x, PNew(1,:), solNew.x);
    parameters.PY = interp1(solOld.x, PNew(2,:), solNew.x);
    
%     parameters.K = K.*gammaOld;
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (trapz(solOld.x, gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
        
        Sols = Sols(1:(i - 1));
        gammaSols = gammaSols(1:(i - 1));
        times = times(1:(i - 1));
        
        break
    end
    
    solOld = solNew;

    Sols{i} = [solOld.x; solOld.y];
    gammaSols{i} = [L.*solNew.x; gammaOld];
                        
end

toc

% Save the solutions

outputDirectory = '../../Solutions/LinearViscoelasticFoundation/';    
save([outputDirectory, 'sols_Eb_0p75_nu_0p02_k_0p02_L0_0p125_sigma_0p1L_area_1_mu_0_inext.mat'], 'Sols') % Solutions
save([outputDirectory, 'gamma_Eb_0p75_nu_0p02_k_0p02_L0_0p125_sigma_0p1L_area_1_mu_0_inext.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'times_Eb_0p75_nu_0p02_k_0p02_L0_0p125_sigma_0p1L_area_1_mu_0_inext.mat'], 'times') % Times
