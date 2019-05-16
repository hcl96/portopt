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

for rhoBS = -1:1:1
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

%% Chapter 3 Demo: 2 Asset Efficient Frontier
clear;clc;
stdS = 12.2;
stdB = 5.5;
returnS = 10.3;
returnB = 6.2;
rhoBS = 0.34;

% create plot of ER and varR by varying Xb
dx = 0.01
xb = 0:dx:1; returnALL = []; varALL = [];
% including short sell: any Xb
for i = 1:1:size(xb,2)
    returnALL = [returnALL xb(i)*returnB + (1 - xb(i))*returnS];
    varALL = [varALL sqrt((xb(i)*stdB)^2 + ((1 - xb(i))*stdS)^2 + 2*xb(i)*(1 - xb(i))*stdB*stdS*rhoBS)];
end
figure(1)
plot(varALL, returnALL); hold on;
title("Risk-Return Plot: Perfect Correlation")
xlabel("Standard Deviation")
ylabel("Expected Return")

%% Chapter 4 Demo: Multi-Asset Efficient Frontier (IN PROGRESS)
clear;clc;
% returns/std for SP500,bonds,CAD,JPN,Emerging,Pacific,EUR,Small Stocks
returnList = [14 6.5 11 14 16 18 12 17];
stdList = [18.5 5 16 23 30 26 20 24];
[returnList;stdList]
figure(1)
plot(stdList, returnList, '*'); hold on;
CorrMat = [1 0.45 0.7 0.2 0.64 0.3 0.61 0.79;
    0 1 0.27 -0.01 0.41 0.01 0.13 0.28;
    0 0 1 0.14 0.51 0.29 0.48 0.59;
    0 0 0 1 0.25 0.73 0.56 0.13;
    0 0 0 0 1 0.28 0.61 0.75;
    0 0 0 0 0 1 0.54 0.16;
    0 0 0 0 0 0 1 0.44;
    0 0 0 0 0 0 0 1];
CorrMat = (CorrMat+CorrMat') - eye(size(CorrMat,1)).*diag(CorrMat); % Make symmetric matrix
CovMatPre = diag(stdList.^2-1)+ones(max(size(stdList)));
CovMat = CovMatPre.*CorrMat
% solve for minimal variance given target return
VarALL = []; ReturnALL = [];
dr = 0.01;
r = 0:dr:25;
C = CovMat; Rbar = returnList'; e = ones(8,1);
for i = 1:1:max(size(r))
    iterConst = inv([e'*inv(C)*Rbar e'*inv(C)*e; Rbar'*inv(C)*Rbar Rbar'*inv(C)*e])*[1;r(i)];
    alpha = iterConst(1); beta = iterConst(2);
    x = alpha*inv(C)*Rbar + beta*inv(C)*e;
    VarALL = [VarALL x'*C*x];
    ReturnALL = [ReturnALL r(i)];
end
plot(sqrt(VarALL), ReturnALL)
title("Risk-Return Plot: Multi-Asset")
xlabel("Standard Deviation")
ylabel("Expected Return")















