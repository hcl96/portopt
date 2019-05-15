%% Problem 1.8
clear;clc;
prob = [0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1];
returns = [-3, -2, -1, 0, 1, 2, 10];
prob_return = [prob;returns]'
meanReturn = prob*returns'
varReturn = 0;
for i = 1:size(returns, 2)
    varReturn = varReturn + prob(i)*(returns(i) - meanReturn)^2;
end
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
% var_returnA = xA*Cov*xA'
% var_returnB = xB*Cov*xB'
pre = diag([1/sqrt(6),1/sqrt(5),1/sqrt(4)])
post = pre
pre*Cov*post %correlation of retrurns given by the off-diagonal terms

%% Chapter 2 Demo: 2 Asset Portfolio
clear;clc;
stdS = 0.15;
stdB = 0.1;
returnS = 0.06;
returnB = 0.05;
%%% ------- Correlation Parameter: Change to tune plot -------- %%%
rhoBS = 0;
%%% ----------------------------------------------------------- %%%
for rhoBS = -1:0.1:1
    % create plot of ER and varR by varying Xb
    dx = 0.1
    xb = -10:dx:10; returnALL = []; varALL = [];
    % including short sell: any Xb
    for i = 1:1:size(xb,2)
        returnALL = [returnALL xb(i)*returnB + (1 - xb(i))*returnS];
        varALL = [varALL sqrt((xb(i)*stdB)^2 + ((1 - xb(i))*stdS)^2 + 2*xb(i)*(1 - xb(i))*stdB*stdS*rhoBS)];
    end
    figure(1)
    plot(varALL, returnALL); hold on;
    % no short sell: Xb+Xs=1 strictly enforced
    dx = 0.1
    xb = 0:dx:1; returnALL = []; varALL = [];
    for i = 1:1:size(xb,2)
        returnALL = [returnALL xb(i)*returnB + (1 - xb(i))*returnS];
        varALL = [varALL sqrt((xb(i)*stdB)^2 + ((1 - xb(i))*stdS)^2 + 2*xb(i)*(1 - xb(i))*stdB*stdS*rhoBS)];
    end
    plot(varALL, returnALL);
    title("Risk-Return Plot: Perfect Correlation")
    xlabel("Standard Deviation")
    ylabel("Expected Return")
    pause(1)
end




