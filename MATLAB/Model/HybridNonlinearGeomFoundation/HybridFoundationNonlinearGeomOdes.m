function Odes = HybridFoundationNonlinearGeomOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
nu1 = parameters.nu1;
nu3 = parameters.nu3;
dt = parameters.dt;
gamma = parameters.gamma;
Es = parameters.Es;
Eb = parameters.Eb;
ext = parameters.ext;

r1Old = MOld.y(2,:);
r3Old = MOld.y(3,:);
P1Old = parameters.P1;
P3Old = parameters.P3;

    function dMdS = HybridFoundationEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        r1 = M(2,:);
        r3 = M(3,:);
        n1 = M(4,:);
        n3 = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Es) > 1)
            Es = interp1(MOld.x, Es, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(MOld.x, Eb, x);
        end
        
        r1Old = interp1(MOld.x, r1Old, x);
        r3Old = interp1(MOld.x, r3Old, x);
        P1Old = interp1(MOld.x, P1Old, x);
        P3Old = interp1(MOld.x, P3Old, x);
        
        % If the model is extensible, set alpha to scale with axial tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + n3./Es;
        else
            alpha = 1;
        end
        
        dSdS = L.*ones(1, length(S));
        dr1dS = -L.*r3.*gamma.*m./Eb;
        dr3dS = L.*(gamma.*alpha + r1.*gamma.*m./Eb);
%         dn1dS = L*(K.*alpha.*gamma.*(r1 - (P1Old + dt*nu1.*(r1Old - P1Old))) - n3.*gamma.*m./Eb);
        dn1dS = L*(K.*alpha.*gamma.*(r1 - r1Old + P1Old.*(1 - dt*nu1)) - n3.*gamma.*m./Eb);
        dn3dS = L*(K.*alpha.*gamma.*(r3 - r3Old + P3Old.*(1 - dt*nu3)) + n1.*gamma.*m./Eb);
        dthetadS = L.*gamma.*m./Eb;
        dmdS = -L.*gamma.*alpha.*n1;

        dMdS = [dSdS; dr1dS; dr3dS; dn1dS; dn3dS; dthetadS; dmdS];
        
    end

Odes = HybridFoundationEqns(x, M);

end
