%% Problem 1.8
clear;clc;
prob = [0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1];
returns = [-3, -2, -1, 0, 1, 2, 10];
prob_return = [prob;returns]'
meanReturn = prob*returns'
varReturn = 0;
for i = 1:size(returns, 2)
    varReturn = varReturn + (prob(i)*returns(i)^2);
end
varReturn = varReturn - meanReturn
stdReturn = sqrt(varReturn)
%% Problem 1.8
clear;clc;
exp_return = [10, 12, 14];
Cov = [6,1,0;1,5,1;0,1,4];
xA = [0.5, 0.5, 0];
xB = [0, 0.5, 0.5];
% find expected returns of A and B
exp_returnA = xA*exp_return'
exp_returnB = xB*exp_return'
var_returnA = xA*Cov*xA'
var_returnB = xB*Cov*xB'
