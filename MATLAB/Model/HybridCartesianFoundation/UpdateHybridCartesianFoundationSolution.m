function [SolNew, gammaNew, PNew] = UpdateHybridCartesianFoundationSolution(solMesh, solOld, W, parameters, options)

% Get the relevant parameters to update growth and the foundation in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);

n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

n3s = parameters.n3s;
sigma = parameters.sigma;
mu = parameters.mu;
nu1 = parameters.nu1;
nu3 = parameters.nu3;
dt = parameters.dt;

% Spring stresses
PXOld = parameters.PX;
PYOld = parameters.PY;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));

% Define the ODEs
Odes = @(x, M) HybridFoundationCartesianOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

XNew = SolNew.y(2,:);
YNew = SolNew.y(3,:);
thetaNew = SolNew.y(6,:);

PXNew = 0.5*XNew.*(1 + cos(2*thetaNew)) + 0.5*YNew.*sin(2*thetaNew) + (PXOld - 0.5*XOld.*(1 + cos(2*thetaOld)) - 0.5*YOld.*sin(2*thetaOld)) ...
    + (0.5*XOld.*cos(2*thetaOld) - 0.5*YOld.*cos(2*thetaOld) - PYOld).*(thetaNew - thetaOld) ...
    + dt*(0.5*nu1*(XOld - PXOld).*(1 - cos(2*thetaOld)) - 0.5*nu1*(YOld - PYOld).*sin(2*thetaOld) ...
    - 0.5*nu3*PXOld.*(1 + cos(2*thetaOld)) - 0.5*nu3*PYOld.*sin(2*thetaOld));
PYNew = 0.5*YNew.*(1 - cos(2*thetaNew)) + 0.5*XNew.*sin(2*thetaNew) + (PYOld - 0.5*YOld.*(1 - cos(2*thetaOld)) - 0.5*XOld.*sin(2*thetaOld)) ...
    - (0.5*XOld.*sin(2*thetaOld) - 0.5*YOld.*cos(2*thetaOld) - PXOld).*(thetaNew - thetaOld) ...
    + dt*(0.5*nu1*(YOld - PYOld).*(1 + cos(2*thetaOld)) - 0.5*nu1*(XOld - PXOld).*sin(2*thetaOld) ...
    - 0.5*nu3*PYOld.*(1 + cos(2*thetaOld)) - 0.5*nu3*PXOld.*sin(2*thetaOld));

PNew = [PXNew; PYNew];
