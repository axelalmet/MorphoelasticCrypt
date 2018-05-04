function [SolNew, gammaNew, PNew] = UpdateNonlinearGeomLinearViscoelasticSolution(solMesh, solOld, W, parameters, options)

% Get the relevant parameters to update growth and the foundation in time
r1Old = solOld.y(2,:);
r3Old = solOld.y(3,:);
n3Old = solOld.y(5,:);

n3s = parameters.n3s;
sigma = parameters.sigma;
mu = parameters.mu;
nu1 = parameters.nu1;
nu3 = parameters.nu3;
dt = parameters.dt;

% Spring stresses
P1Old = parameters.P1;
P3Old = parameters.P3;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));

% Define the ODEs
Odes = @(x, M) HybridFoundationNonlinearGeomOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonlinearGeomBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

r1New = SolNew.y(2,:);
r3New = SolNew.y(3,:);

P1New = r1New - r1Old - (1 - dt*nu1).*P1Old;
P3New = r3New - r3Old - (1 - dt*nu3).*P3Old;

PNew = [P1New; P3New];
