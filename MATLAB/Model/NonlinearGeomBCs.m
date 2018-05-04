function Residuals = NonlinearGeomBCs(Ml, Mr, parameters)

% Here:
% M(1) = S, M(2) = r1, M(3) = r3, M(4) = n1, M(5) = n3, 
% M(6) = theta, M(7) = theta;

L = parameters.L;

Residuals = [Mr(1) - L, Ml(2), Mr(2), Ml(3), Mr(3) - L, Ml(6), Mr(6)];
