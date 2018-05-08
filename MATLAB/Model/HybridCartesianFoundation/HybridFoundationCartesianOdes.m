function Odes = HybridFoundationCartesianOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
nu1 = parameters.nu1;
nu3 = parameters.nu3;
dt = parameters.dt;
gamma = parameters.gamma;
Es = parameters.Es;
Eb = parameters.Eb;
ext = parameters.ext;

XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
PXOld = parameters.PX;
PYOld = parameters.PX;

    function dMdS = HybridFoundationEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
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
        
        XOld = interp1(MOld.x, XOld, x);
        YOld = interp1(MOld.x, YOld, x);
        PXOld = interp1(MOld.x, PXOld, x);
        PYOld = interp1(MOld.x, PYOld, x);
        
        % If the model is extensible, set alpha to scale with axial tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        dSdS = L.*ones(1, length(S));
        dxdS = -L.*alpha.*gamma.*cos(theta);
        dydS = L.*gamma.*alpha.*sin(theta);
        dFdS = L*(K.*alpha.*gamma.*(r1 - r1Old + P1Old.*(1 - dt*nu1)) - n3.*gamma.*m./Eb);
        dGdS = L*(K.*alpha.*gamma.*(r3 - r3Old + P3Old.*(1 - dt*nu3)) + n1.*gamma.*m./Eb);
        dthetadS = L.*gamma.*m./Eb;
        dmdS = L.*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = HybridFoundationEqns(x, M);

end
