function Odes = KellerModelContactOdes(x, M, region, modelParams)

P = modelParams.P;
L = modelParams.L;

    function dMdS = ContactEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        x = M(2,:);
        y = M(3,:);
        Q = M(4,:);
        N = M(5,:);
        theta = M(6,:);
        k = M(7,:);
        sc = M(8,:);
        fc = M(9,:);
        
        dsCdS = zeros(1, length(Q));
        dfCdS = zeros(1, length(Q));
        
        switch region
            case 1 % ODEs for [0, sc]
                
                dSdS = (sc).*ones(1, length(S));
                dxdS = (sc).*cos(theta);
                dydS = (sc).*sin(theta);
                dQdS = -(sc).*k.*N;
                dNdS = (sc).*(k.*Q + P);
                dthetadS = (sc).*k;
                dkdS = (sc).*N;
                
            case 2 % ODEs for [sc, L]
                
                dSdS = (L - sc).*ones(1, length(S));
                dxdS = (L - sc).*cos(theta);
                dydS = (L - sc).*sin(theta);
                dQdS = -(L - sc).*(k.*N);
                dNdS = (L - sc).*(k.*Q + P);
                dthetadS = (L - sc).*k;
                dkdS = (L - sc).*N;
                
            otherwise
                error('MATLAB:contactOdes:BadRegionIndex','Incorrect region index: %d',region);
        end
                
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS; dsCdS; dfCdS];
    end

Odes = ContactEqns(x, M, region);

end
