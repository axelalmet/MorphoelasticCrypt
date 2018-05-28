function [SolNew, gammaNew, PNew] = UpdateLinearKelvinVoigtSolution(solMesh, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);

DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

n3s = parameters.n3s;
sigma = parameters.sigma;
mu = parameters.mu;
nu = parameters.nu;
dt = parameters.dt;

% Spring stresses
y0 = parameters.y0;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));

% Define the ODEs
Odes = @(x, M) LinearKelvinVoigtFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

SNew = SolNew.y(1,:);
XNew = SolNew.y(2,:);
YNew = SolNew.y(3,:);
DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

PNew = DeltaOld - y0  + (dt*nu).^(-1).*(DeltaNew - DeltaOld);
