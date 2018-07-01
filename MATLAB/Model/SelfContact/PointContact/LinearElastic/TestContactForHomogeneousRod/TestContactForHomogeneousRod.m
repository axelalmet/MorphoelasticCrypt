function TestContactForHomogeneousRod
% Set the parameters
h = 0.015; % Thickness of the rod cross section
L0 = 0.125; % Dimensional length of the rod`
L = 2*sqrt(3)*L0/h; %Dimensionless length
y0 = 0;
K = 0.01; % Dimensionless foundation stiffness
Es = 1; % Stretching stiffness
Eb = 1; % Bending stiffness
dt = 5*1e-2;
w0 = 0;

% Get the initial solution from AUTO
solData = load('../../../../../../Data/planarmorphorodsk0p01L29halfint_sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

solMesh = solFromData.x;
 
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
parameters.y0 = y0;
parameters.Es = Es; % Stretch stiffness
parameters.Eb = Eb; % Bending stiffness
parameters.ext = 0; % Extensibility
parameters.dt = dt; % Time step
parameters.w0 = w0;

%% Solve the BVP until contact

parameters.gamma = 1 + dt;

tic
while (parameters.gamma < 2.1)

% Define the ODEs and BCs
DerivFun = @(x, M) LinearElasticFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) HalfIntervalBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);

parameters.gamma = (parameters.gamma) + dt;

solFromData.x = solMesh;
solFromData.y = deval(numSol, solMesh);

end

toc

%% Construct the initial guess for bvp4c 

% Initialise the guesses for the contact points
possibleContactPoints = find(abs(solFromData.y(2,:) - 0.5*L) < 0.125);
contIndexOne = possibleContactPoints(1);
contIndexTwo = possibleContactPoints(end);
sCGuess = solFromData.y(1, contIndexOne);
fCGuess = 0;

initSol.x = [solFromData.x(1:contIndexOne), solFromData.x((contIndexOne):end)];

initSol.y = [solFromData.y(:, 1:contIndexOne), solFromData.y(:, contIndexOne:end); ...
    sCGuess.*ones(1, length(initSol.x)); fCGuess.*ones(1, length(initSol.x))];

initSol.x = [linspace(0, 1, contIndexOne), linspace(1, 2, length(initSol.x) - contIndexOne)];
initSol.y = initSol.y(:, 1:(end - 1)/2);

initSol.y(1:2, end) = [0.5*L; 0.5*L];
% 
solMesh = initSol.x;

%% Solve 

parameters.gamma = 2.01 + dt;

% Define the ODEs and BCs
DerivFun = @(x, M, region) LinearElasticFoundationContactOdes(x, M, region, initSol, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) SelfPointContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, initSol, solOptions);

firstContactSol = numSol;

toc

%% Grow the rod with contact

while ( (parameters.gamma < 7.0 )&&(numSol.y(9,1) > 0) )

% Define the ODEs and BCs
DerivFun = @(x, M, region) LinearElasticFoundationContactOdes(x, M, region, numSol, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) SelfPointContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, numSol, solOptions);

parameters.gamma = (parameters.gamma) + dt;

numSol.x = solMesh;
numSol.y = deval(newSol, solMesh);

end

toc