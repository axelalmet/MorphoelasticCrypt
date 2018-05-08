using DifferentialEquations
using Interpolations

# Import initial guess
initSolData = readdlm("../Data/planarmorphorodsk0p02L29_sol_1")
solMesh = initSolData[:,1]'
initS = initSolData[:,2]'
initX = initSolData[:,3]'
initY = initSolData[:,4]'
initF = initSolData[:,5]'
initG = initSolData[:,6]'
initTheta = initSolData[:,7]'
initM = initSolData[:,8]'

initR1 = -initX.*sin(initTheta) + initY.*cos(initTheta)
initR3 = initX.*cos(initTheta) + initY.*sin(initTheta)
initN1 = -initF.*sin(initTheta) + initG.*cos(initTheta)
initN3 = initF.*cos(initTheta) + initG.*sin(initTheta)

# Define the parameters
kf = 0.16 # Dimensional foundational stiffness
h = 0.015 # Thickness of the rod cross section
w = 0.01 # Width of the rod cross section
L0 = 0.125 # Dimensional length of the rod
L = 2.0*sqrt(3.0)*L0/h  # Dimensionless length
K = kf*h/(12.0*w) # Dimensionless foundation stiffness
n3s = 0.0 # Target axial tension
Es = 1.0 # Stretching stiffness
Eb = 1.0 # Bending stiffness
dt = 1e-4 # Time step

sigma = 0.1*L # Width of Wnt gaussian
W(x) = exp(-(L*(x - 0.5)/sigma).^2.0) # Define Wnt function

# Define the growth sensitivities
eta = quadgk(W, 0.0, 1.0)
mu = 0.0

# Define the viscoelastic foundation parameters
nu1 = K; # Normal direction
nu3 = K; # Tangent direction

ext = 0.0 # Set extensibility assumption

initSol = [initS; initR1; initR3; initN1; initN3; initTheta; initM]

function Odes!(du, u, p, t)

    S = u[1,:]
    r1 = u[2,:]
    r3 = u[3,:]
    n1 = u[4,:]
    n3 = u[5,:]
    theta = u[6,:]
    m = u[7,:]

    if [ext == 1]
        alpha = 1.0 + n3.*Es.^(-1)
    else()
        alpha = 1.0
    end

    dS = L*ones(S)
    dr1 = -L.*r3*gamma.*m.*Eb.^(-1)
    dr3 = L*(gamma.*alpha + r1.*gamma.*m.*Eb.^(1)
    dn1 = L*(K.*alpha.*gamma.*(r1 - (P1Old + dt*nu1.*(r1Old - P1Old))) - n3.*gamma.*m.*Eb.^(-1))
    dn3 = L*(K.*alpha.*gamma.*(r3 - r3Old + P3Old.*(1 - dt*nu3)) + n1.*gamma.*m.*Eb.^(-1))
    dtheta = L.*gamma.*m.*Eb.^(-1)
    dm = -L.*gamma.*alpha.*n1

    du = [dS; dr1; dr3; dn1; dn3; dtheta; dm]
end

function Bcs!(residual, u, p, t)

    residual[1] = u[end][1] - L # S(L) = L
    residual[2] = u[1][2] # r1(0) = 0
    residual[3] = u[end[2] # r1(L) = 0
    residual[4] = u[1][3] # r3(0) = 0
    residual[5] = u[end][3] - L # r3(L) = L
    residual[6] = u[1][4] # theta(0) = 0
    residual[7] = u[end][4] #theta(L) = 0

end