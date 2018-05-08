function Residuals = NonUniformGrowthBCs(Ml, Mr, parameters)

L = parameters.L;

Residuals = [Ml(1), Ml(2), Mr(2) - L, Ml(3), Mr(3), Ml(6), Mr(6)];
