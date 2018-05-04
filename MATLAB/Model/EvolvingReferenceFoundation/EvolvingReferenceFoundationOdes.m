function Odes = EvolvingReferenceFoundationOdes(x, M, MOld, parameters)

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

    function dMdS = EvolvingReferenceEqns(x, M)
        
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
            gamma = interp1(MOld.x, gamma, x).*(length(gamma) > 1);
        end
        
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
        
        if (length(r1Old) > 1)
            r1Old = interp1(MOld.x, r1Old, x);
        end
        
        if (length(r3Old) > 1)
            r3Old = interp1(MOld.x, r3Old, x);
        end
        
        if (length(P1Old) > 1)
            P1Old = interp1(MOld.x, P1Old, x);
        end
        
        if (length(P3Old) > 1)
            P3Old = interp1(MOld.x, P3Old, x);
        end
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + n3./Es;
        else
            alpha = 1;
        end
        
        dSdS = L.*ones(1, length(S));
        dr1dS = zeros(1, length(S));
        dr3dS = L.*gamma.*alpha;
        dn1dS = L*K.*alpha.*gamma.*(P1Old + dt*nu1.*(r1Old - P1Old));
        dn3dS = L*K.*alpha.*gamma.*(P3Old + dt*nu3.*(r3Old - P3Old));
        dthetadS = L.*gamma.*m./Eb;
        dmdS = -L.*gamma.*alpha.*n1;

        dMdS = [dSdS; dr1dS; dr3dS; dn1dS; dn3dS; dthetadS; dmdS];
        
    end

Odes = EvolvingReferenceEqns(x, M);

end
