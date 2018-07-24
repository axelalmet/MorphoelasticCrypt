function Residuals = KellerModelContactRegionBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact along a region

L = modelParams.L;

Residuals = [Ml(1, 1), ...
             Ml(2, 1), ... x(0) = 0
             Ml(3, 1), ... y(0) = 0
             Ml(5, 1), ...
             Ml(6, 1), ... theta(0) = 0
             Ml(1, 2) - Mr(1, 1), ... S(s1+) = S(s1-)
             Mr(2, 1), ... x(s1) = 0
             Ml(3, 2) - Mr(3, 1), ... y(s1+) = y(s1-)
             Mr(4, 1), ...
             Ml(4, 2) - (Mr(4, 1) - Mr(9, 1)), ... Q(s1+) = Q(s1-) - fc
             Mr(5, 1), ...
             Ml(5, 2) - Mr(5, 1), ... N(s1+) = N(s1-)
             Mr(6, 1) - 0.5*pi, ... theta(s1-) = pi/2
             Mr(7, 1), ... k(s1) = 0
             Ml(9, 2) - Mr(9, 1), ... % fC(s1+) = fC(s1-)
             Mr(1, 2) - L];
         
end