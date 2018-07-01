function [solMeshNew, SolNew, gammaNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);
mOld = solOld.y(7,:);

n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

n3s = parameters.n3s;
sigma = parameters.sigma;
mu = parameters.mu;
nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
b1 = parameters.b1;
Eb = parameters.Eb;

% Spring stresses
POld = parameters.P;
uHatOld = parameters.uHat;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));

% Set new bending stiffness
parameters.Eb = 1 - b1.*W(solOld.x, sigma./gammaOld);

% Define the ODEs
Odes = @(x, M) SimplifiedLinearMaxwellFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

distanceMatrix = abs(Sol.y(2,2:end) - Sol.y(2,1:(end - 1)));

% Re-mesh the new solution if the inf-norm is too large
if (max(distanceMatrix(:)) > 0.25) 
    
    solMeshNew = linspace(0, 1, 2*length(solMeshOld));   
    
    % For sanity
    solMeshNew(1) = 0;
    solMeshNew(end) = 1;
    
else
    
    solMeshNew = solMeshOld;
    
end

SolNew.x = solMeshNew;
SolNew.y = deval(Sol, solMeshNew);

YNew = interp1(solMeshNew, SolNew.y(3,:), solMeshOld);

PNew = YNew - YOld + (1 - dt*nu).*POld; 

uHatNew = uHatOld + etaK*dt.*mOld./Eb;