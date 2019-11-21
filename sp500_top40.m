%% Find Mean Return
% Data taken from https://finance.yahoo.com
% 5/1/2009 ~ 5/1/2019, Monthly Close
% Return Time Horizon = 1 year
% Stock without Info: BRKB
% Stock that has insufficient Info: FB, PYPL
clear;clc;
stockList = {'MSFT','AAPL','AMZN','JNJ','JPM','GOOG','GOOGL','XOM'...
    ,'V','PG','BAC','CSCO','VZ','UNH','DIS','PFE','T','MA','CVX','HD',...
    'MRK','INTC','KO','CMCSA','WFC','BA','PEP','NFLX','MCD','C','WMT',...
    'ABT','ADBE','ORCL','MDT','PM','UNP'};
stock_close = [];
for i = 1:max(size(stockList))
    disp(stockList{i})
    raw = readtable(strcat(stockList{i},'.csv'));
    if size(raw,1) == 121
        stock_close = [stock_close table2array(raw(:,5))];
    else
        disp('size is small!')
    end
    
end
figure(1)
plot(stock_close(:,:));
title('Top 10 S&P500 Stocks by Weight ');
xlabel('Time (Day)'); ylabel('Price ($)');
legend(stockList)
%% Find Expected Return
yearly_ret = []; monthly_ret = []; retrec = [];
for i = 1:size((stock_close),2)
    ret = [];
    for j = 1:size((stock_close),1)-1
        % find daily return rate, then find mean daily return
        ret = [ret (stock_close(j+1,i) - stock_close(j,i))/stock_close(j,2)];
    end
    % yearly return rate is daily return rate compounded over a year
    yearly_ret = [yearly_ret 100*((1+mean(ret))^(12)-1)];
    retrec = [retrec; ret*100];
    monthly_ret = [monthly_ret mean(ret)*100];
end
%% Find Covariance Matrix
CovMat = zeros(size((stock_close),2),size((stock_close),2));
for i = 1:size((stock_close),2)
    for j = 1:size((stock_close),2)
        for k = 1:size((stock_close),1)-1
            CovMat(i,j) = CovMat(i,j) + (retrec(i,k)-monthly_ret(i))*(retrec(j,k)-monthly_ret(j));
        end
        CovMat(i,j) = CovMat(i,j)/(size((stock_close),1)-1);
    end
end

%%
figure(2)
title('Mean Daily Return for MSFT')
xlabel('Time (Day)'); ylabel('Price ($)');
plot(retrec(1,:)); hold on;
plot(mean(retrec(1,:))*ones(size(stock_close,1),1)); hold off;

figure(3) % plotting mean and variance for daily return
n_hist = 10; % number of histogram bins
subplot(2,2,1)
histogram(retrec(1,:), n_hist)
title('Expected Return Distribution- MFST')
subplot(2,2,2)
histogram(retrec(2,:), n_hist)
title('Expected Return Distribution- AAPL')
subplot(2,2,3)
histogram(retrec(3,:), n_hist)
title('Expected Return Distribution- FB')
subplot(2,2,4)
histogram(retrec(4,:), n_hist)
title('Expected Return Distribution- JPM')
%%
returnList = monthly_ret;
stdList = sqrt(diag(CovMat));
% solve for minimal variance for ANY LEVERAGE
stdAnyLeverage = []; stdNoLeverage = []; 
ReturnAnyLeverage = []; ReturnNoLeverage = [];

xrecAnyLeverage = []; xrecNoLeverage = []; 
for r = 0:0.1:15
    q = 2*CovMat; % minimize 0.5*x'Qx
    c = [];
    a = [returnList; ones(1,size(returnList,2))];
    blc = [r; 1]; % lower limit is target return
    buc = [r; 1]; % upper limit is target return
    blx = [];
    bux = [];
    [res] = mskqpopt(q,c,a,blc,buc,blx,bux);
    x = res.sol.itr.xx;
    stdAnyLeverage = [stdAnyLeverage sqrt(x'*CovMat*x)];
    ReturnAnyLeverage = [ReturnAnyLeverage r];
    xrecAnyLeverage = [xrecAnyLeverage x]
end
% solve for minimal variance for NO LEVERAGE
x = [];
for r = 0:0.1:15
    q = 2*CovMat; % minimize 0.5*x'Qx
    c = [];
    a = [returnList; ones(1,size(returnList,2))];
    blc = [r; 1]; % lower limit is target return
    buc = [r; 1]; % upper limit is target return
    blx = sparse(size(returnList,2),1);
    bux = [];
    [res] = mskqpopt(q,c,a,blc,buc,blx,bux);
    x = res.sol.itr.xx;
    if mean(x) ~= 0
        stdNoLeverage = [stdNoLeverage sqrt(x'*CovMat*x)];
        ReturnNoLeverage = [ReturnNoLeverage r];
        xrecNoLeverage = [xrecNoLeverage x]
    end
end
%%
figure(1)
plot(stdList, returnList, 'bx'); hold on;
plot(stdAnyLeverage, ReturnAnyLeverage, 'b--')
plot(stdNoLeverage, ReturnNoLeverage, 'r-')
title("Top 40 S&P 500 Stocks Efficient Frontier")
xlabel("Standard Deviation (Risk)")
ylabel("Expected Return")
text(stdList, returnList,stockList,...
    'VerticalAlignment','top','HorizontalAlignment','left')

%% find optimum portfolio
targetReturn = 4; %
k = find(ReturnNoLeverage == targetReturn);
xrecNoLeverage(:, k)


