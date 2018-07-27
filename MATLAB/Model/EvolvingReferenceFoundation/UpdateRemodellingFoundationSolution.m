function [solMeshNew, SolNew, gammaNew, EbNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolution(solMeshOld, solOld, W, parameters, options)
% Get the relevant parameters to update growth in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

mu1 = parameters.mu1;
mu2 = parameters.mu2;
width = parameters.width;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;
nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
b1 = parameters.b1;
Eb = parameters.Eb;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;
uHatOld = parameters.uHat;

% Define new gamma
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*W(solOld.y(1,:), sigma));
gammaNew = gammaOld.*(1 + dt.*W(currentArcLength, mu1, mu2, width));

% Set new bending stiffness
% EbNew = 1 - b1.*W(currentArcLength, mu1, mu2, width);
EbNew = 1 - b1.*exp(-((currentArcLength - 0.5*(mu1 + mu2))./sigma).^2);

parameters.Eb = EbNew;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solOld, parameters);

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

PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

uHatNew = uHatOld + etaK*dt.*mOld./Eb;