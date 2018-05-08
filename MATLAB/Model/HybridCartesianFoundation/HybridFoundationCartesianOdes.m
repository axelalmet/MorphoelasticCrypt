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
thetaOld = MOld.y(6,:);
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
        thetaOld = interp1(MOld.x, thetaOld, x);
        PXOld = interp1(MOld.x, PXOld, x);
        PYOld = interp1(MOld.x, PYOld, x);
        
        % If the model is extensible, set alpha to scale with axial tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        % Update the spring stresses
        PX = 0.5*X.*(1 + cos(2*theta)) + 0.5*Y.*sin(2*theta) + (PXOld - 0.5*XOld.*(1 + cos(2*thetaOld)) - 0.5*YOld.*sin(2*thetaOld)) ...
            + (0.5*XOld.*cos(2*thetaOld) - 0.5*YOld.*cos(2*thetaOld) - PYOld).*(theta - thetaOld) ...
            + dt*(0.5*nu1*(XOld - PXOld).*(1 - cos(2*thetaOld)) - 0.5*nu1*(YOld - PYOld).*sin(2*thetaOld) ...
                - 0.5*nu3*PXOld.*(1 + cos(2*thetaOld)) - 0.5*nu3*PYOld.*sin(2*thetaOld));
        PY = 0.5*Y.*(1 - cos(2*theta)) + 0.5*X.*sin(2*theta) + (PYOld - 0.5*YOld.*(1 - cos(2*thetaOld)) - 0.5*XOld.*sin(2*thetaOld)) ...
            - (0.5*XOld.*sin(2*thetaOld) - 0.5*YOld.*cos(2*thetaOld) - PXOld).*(theta - thetaOld) ...
            + dt*(0.5*nu1*(YOld - PYOld).*(1 + cos(2*thetaOld)) - 0.5*nu1*(XOld - PXOld).*sin(2*thetaOld) ...
                - 0.5*nu3*PYOld.*(1 + cos(2*thetaOld)) - 0.5*nu3*PXOld.*sin(2*thetaOld));
        
        % Model derivatives
        dSdS = L.*ones(1, length(S));
        dxdS = L.*alpha.*gamma.*cos(theta);
        dydS = L.*gamma.*alpha.*sin(theta);
        dFdS = L*(K.*alpha.*gamma.*(PX.*cos(2*theta) + PY.*sin(2*theta) ...
                + 0.5*X.*(1 - cos(2*theta)) - 0.5.*Y.*sin(2*theta)) ...
                - G.*gamma.*m./Eb);
        dGdS = L*(K.*alpha.*gamma.*(PX.*sin(2*theta) - PY.*cos(2*theta) ...
                - 0.5*X.*sin(2*theta) + 0.5*Y.*(1 + cos(2*theta))) + F.*gamma.*m./Eb);
        dthetadS = L.*gamma.*m./Eb;
        dmdS = L.*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = HybridFoundationEqns(x, M);

end
