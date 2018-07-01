function Residuals = SelfPointContactBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact at a single point.
% Specifically, we split the problem into three intervals and introduce
% matching conditions. This results in 21 boundary conditions.

% Obtain the relevant model parameters
w0 = modelParams.w0;
y0 = modelParams.y0;
L = modelParams.L;

Residuals = [Ml(1, 1), ... % S1(0) = 0
             Ml(2, 1), ... % x1(0) = 0
             Ml(3, 1) - y0, ... % y1(0) = y0
             Ml(6, 1), ... % theta1(0) = 0
             Mr(2, 1) - (L - w0), ... % x1(S1) = L/2 - w; contact point
             Ml(2, 2) - (L - w0), ... % x2(S1) = L/2 - w; contact point
             Mr(6, 1) - 0.5*pi, ... % theta1(S1) = pi/2; contact point)
             Ml(6, 2) - 0.5*pi, ... % theta2(s1) = pi/2
             Ml(4, 2) - (Mr(4, 1) + Mr(9, 1)), ... % F2(S1) = F1(S1) + fC; jump condition
             Ml(5, 2) - Mr(5, 1), ... % G2(S1) = G1(S2); continuity
             Ml(1, 2) - Mr(1, 1), ... % G2(S1) = G1(S2); continuity
             Ml(3, 2) - Mr(3, 1), ... % y3(S2) = y2(S2)
             Ml(7, 2) - Mr(7, 1), ... % m(S2) = m(S2)
             Ml(9, 2) - Mr(9, 1), ...  % enforce symmetric contact
             Mr(1, 2) - L, ... % S3(L) = L
             Mr(2, 2) - L, ... % x3(L) = L
             Mr(5, 2), ...
             Mr(6, 2)]; % theta3(L) = 0
         
end