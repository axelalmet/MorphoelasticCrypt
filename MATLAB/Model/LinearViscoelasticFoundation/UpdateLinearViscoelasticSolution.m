function [SolNew, gammaNew, PNew] = UpdateLinearViscoelasticSolution(solMesh, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

n3s = parameters.n3s;
sigma = parameters.sigma;
mu = parameters.mu;
eta = parameters.eta;
nu = parameters.nu;
dt = parameters.dt;

% Spring stresses
PXOld = parameters.PX;
PYOld = parameters.PY;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));

% Define the ODEs
Odes = @(x, M) LinearViscoelasticFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

XNew = SolNew.y(2,:);
YNew = SolNew.y(3,:);

PXNew = XNew - XOld + (1 - dt*nu).*PXOld;
PYNew = YNew - YOld + (1 - dt*nu).*PYOld;

PNew = [PXNew; PYNew];
