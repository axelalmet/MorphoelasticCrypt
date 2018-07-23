function TestKellerModel
%% Set-up
P = 5.245;
L = 0.5*pi;

parameters.P = P;
parameters.L = L;

% First derive the buckled solution
solMesh = 0:5*1e-3:1;

initSol.x = solMesh;
initSol.y = [L.*solMesh;
             cos(L.*solMesh + 0.5*sin(2.*L.*solMesh)); ...
             sin(L.*solMesh + 0.5*sin(2.*L.*solMesh)); ...
             0.*solMesh - 0.5*2*cos(2.*L.*solMesh); ...
             0.*solMesh - 0.5*4*sin(2.*L.*solMesh); ...
             L.*solMesh + 0.5*sin(2.*L.*solMesh); ...
             ones(1, length(solMesh)) + 0.5*2*cos(L.*2.*solMesh)];

% Define the ODEs and BCs
DerivFun = @(x, M) KellerModelPreContactOdes(x, M, parameters);
    
% Set the boundary conditions 
BcFun = @(Ml, Mr) KellerModelPreContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, initSol, solOptions);

numSol.x = solMesh;
numSol.y = deval(newSol, solMesh);

toc

%% Do a simple continuation step to contact

contIndex = 180;
% sCGuess = trapz(L.*solMesh(1:contIndex), sqrt(xSol(1:contIndex).^2 + ySol(1:contIndex).^2));
sCGuess = numSol.y(1, 180);

oldSol = [numSol.y(:, 1:contIndex), numSol.y(:, contIndex:end)];
solOld.x = [linspace(0, 1, contIndex), linspace(1, 2, length(numSol.y(1,:)) - contIndex + 1)];
solMesh = solOld.x;
solOld.y = [oldSol; sCGuess.*ones(1, length(solOld.x)); zeros(1, length(solOld.x))];

tic
parameters.P = 5.3;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Define the ODEs and BCs
DerivFun = @(x, M, region) KellerModelContactOdes(x, M, region, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) KellerModelContactBCs(Ml, Mr, parameters);

contactSol = bvp4c(DerivFun, BcFun, solOld, solOptions);

contactSol.x = solMesh;
contactSol.y = deval(contactSol, solMesh);

toc

%% Compare the two solutions
figure(1)
hold on
plot(numSol.y(2,:), numSol.y(3,:))
plot(contactSol.y(2,:), contactSol.y(3,:))

figure(2)
hold on
plot(contactSol.y(2,:), contactSol.y(3,:))

figure(3)
hold on
plot(contactSol.y(1,:), contactSol.y(5,:))

figure(4)
hold on
plot(contactSol.y(1,:), contactSol.y(7,:))

%% Continue the contact solution until we get close to the region contact case

contactSolOld = contactSol;

dP = 0.1;
PMax = 10.5;

tic
while (parameters.P < PMax)
    parameters.P = parameters.P + dP;
    
    % Define the ODEs and BCs
    DerivFun = @(x, M, region) KellerModelContactOdes(x, M, region, parameters);
    
    % Set the boundary conditions
    BcFun = @(Ml, Mr) KellerModelContactBCs(Ml, Mr, parameters);
    
    NewSol = bvp4c(DerivFun, BcFun, contactSolOld, solOptions);
    
    contactSolNew.x = solMesh;
    contactSolNew.y = deval(NewSol, solMesh);
    
    contactSolOld = contactSolNew;
   
    
end
toc

figure(2)
hold on
plot(contactSolOld.y(2,:), contactSolOld.y(3,:))

figure(3)
hold on
plot(contactSolOld.y(1,:), contactSolOld.y(5,:))

figure(4)
hold on
plot(contactSolOld.y(1,:), contactSolOld.y(7,:))