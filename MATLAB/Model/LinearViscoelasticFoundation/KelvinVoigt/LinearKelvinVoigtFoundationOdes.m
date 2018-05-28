function Odes = LinearKelvinVoigtFoundationOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
y0 = parameters.y0;
nu = parameters.nu;
dt = parameters.dt;
gamma = parameters.gamma;
Es = parameters.Es;
Eb = parameters.Eb;
ext = parameters.ext;

SOld = MOld.y(1,:);
XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
POld = parameters.P;

    function dMdS = VicoelasticFoundationEqns(x, M)
        
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
        
        if (length(SOld) > 1)
            SOld = interp1(MOld.x, SOld, x);
        end
        
        if (length(XOld) > 1)
            XOld = interp1(MOld.x, XOld, x);
        end
        
        if (length(YOld) > 1)
            YOld = interp1(MOld.x, YOld, x);
        end
        
        if (length(POld) > 1)
            POld = interp1(MOld.x, POld, x);
        end
       
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        Delta = sqrt((X - S).^2 + (Y).^2);
        DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);
        
        % Update the viscoelastic stress in time
        P = DeltaOld - y0 + (dt*nu)^(-1).*(Delta - DeltaOld); 
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*alpha.*cos(theta);
        dydS = L.*gamma.*alpha.*sin(theta);
        dFdS = L*K.*gamma.*alpha.*P.*(X - S)./Delta;
        dGdS = L*K.*gamma.*alpha.*P.*(Y)./Delta;
%         dFdS = L*K.*P.*(X - S)./Delta;
%         dGdS = L*K.*P.*(Y)./Delta;
        dthetadS = L.*gamma.*m./Eb;
        dmdS = L.*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = VicoelasticFoundationEqns(x, M);

end
