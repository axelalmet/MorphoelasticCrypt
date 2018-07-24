function TestKellerModel
%% Set-up
P = 5.24;
L = 0.5*pi;

parameters.P = P;
parameters.L = L;

% First derive the buckled solution
solMesh = 0:(5*1e-3):1;

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

tic
parameters.P = 5.3;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

solOld = numSol;

% Define the ODEs and BCs
DerivFun = @(x, M) KellerModelContactOdes(x, M, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) KellerModelContactBCs(Ml, Mr, parameters);

contactSol = bvp4c(DerivFun, BcFun, solOld, solOptions);

contactSol.x = solMesh;
contactSol.y = deval(contactSol, solMesh);

toc

%% Compare the two solutions
figure(3)
hold on
plot(numSol.y(2,:), numSol.y(3,:))
plot(contactSol.y(2,:), contactSol.y(3,:))

figure(5)
hold on
plot(contactSol.y(1,:), contactSol.y(7,:))

%% Continue the contact solution until we get close to the region contact case

contactSolOld = contactSol;

dP = 0.02;
PMax = 10.32;

tic
while (parameters.P < PMax)
    parameters.P = parameters.P + dP;
    
    % Define the ODEs and BCs
    DerivFun = @(x, M) KellerModelContactOdes(x, M, parameters);
    
    % Set the boundary conditions
    BcFun = @(Ml, Mr) KellerModelContactBCs(Ml, Mr, parameters);
    
    NewSol = bvp4c(DerivFun, BcFun, contactSolOld, solOptions);
    
    contactSolNew.x = solMesh;
    contactSolNew.y = deval(NewSol, solMesh);
    
    contactSolOld = contactSolNew;
    
end
toc
%%
figure(2)
hold on
plot(contactSolOld.y(2,:), contactSolOld.y(3,:))

figure(3)
hold on
plot(contactSolOld.y(1,:), contactSolOld.y(5,:))

figure(4)
hold on
plot(contactSolOld.y(1,:), contactSolOld.y(7,:))

%% Test for contact along a region

% Construct the new solution guess

% contIndex = 190;
scGuess = L - 0.01;

% oldSol = [contactSolOld.y(:, 1:contIndex), ...
%             [contactSolOld.y(1, contIndex:end); 0.*contactSolOld.y(2,contIndex:end); ...
%              contactSolOld.y(3:5, contIndex:end); 0.5*pi + 0.*contactSolOld.y(2,contIndex:end); ...
%              contactSolOld.y(7, contIndex:end)]];
% solOld.x = [linspace(0, 1, contIndex), ...
%             linspace(1, 2, length(contactSolOld.y(1,:)) - contIndex + 1)];
% solMesh = solOld.x;
% solOld.y = [oldSol; scGuess.*ones(1, length(solOld.x)); zeros(1, length(solOld.x))];
solOld.x = contactSolOld.x;
solOld.y = [contactSolOld.y; scGuess.*ones(1, length(solOld.x))];

tic
parameters.P = 10.4;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Define the ODEs and BCs
% DerivFun = @(x, M, region) KellerModelContactRegionOdes(x, M, region, parameters);
DerivFun = @(x, M) KellerModelSimplifiedContactRegionOdes(x, M, parameters);

% Set the boundary conditions
% BcFun = @(Ml, Mr) KellerModelContactRegionBCs(Ml, Mr, parameters);
BcFun = @(Ml, Mr) KellerModelSimplifiedContactRegionBCs(Ml, Mr);

contactRegionSol = bvp4c(DerivFun, BcFun, solOld, solOptions);

contactRegionSol.x = solMesh;
contactRegionSol.y = deval(contactRegionSol, solMesh);

toc
%% Run a continuation in P to make sure it's doing what it's supposed to!

PMax = 20;
dP = 0.1;

contactRegionSolOld = contactRegionSol;

tic
while (parameters.P < PMax)
    parameters.P = parameters.P + dP;
    
    % Define the ODEs and BCs
    DerivFun = @(x, M) KellerModelSimplifiedContactRegionOdes(x, M, parameters);
    
    % Set the boundary conditions
    BcFun = @(Ml, Mr) KellerModelSimplifiedContactRegionBCs(Ml, Mr);
    
    NewSol = bvp4c(DerivFun, BcFun, contactRegionSolOld, solOptions);
    
    contactRegionSolNew.x = solMesh;
    contactRegionSolNew.y = deval(NewSol, solMesh);
    
    contactRegionSolOld = contactRegionSolNew;
    
    figure(1)
    hold on
    plot(contactRegionSolOld.y(2,:), contactRegionSolOld.y(3,:))
    
    figure(2)
    hold on
    plot(contactRegionSolOld.y(1,:), contactRegionSolOld.y(8,:))
    
end
toc