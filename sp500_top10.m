%% Find Mean Return
% Data taken from https://finance.yahoo.com
% 1/1/2013 ~ 1/1/2019, Daily EOD
% Return Time Horizon = 1 year
clear;clc;
stockList = {'MSFT','AAPL','FB','JPM','AMZN'};
stock_close = [];
for i = 1:max(size(stockList))
    disp(stockList{i})
    raw = readtable(strcat(stockList{i},'.csv'));
    stock_close = [stock_close table2array(raw(:,5))];
end
figure(1)
plot(stock_close(:,1:4));
title('Top 10 S&P500 Stocks by Weight (AMZN Excluded)');
xlabel('Time (Day)'); ylabel('Price ($)');
legend(stockList)

% Find Expected Return
yearly_ret = []; daily_ret = []; retrec = []; varrec = [];
for i = 1:size((stock_close),2)
    ret = [];
    for j = 1:size((stock_close),1)-1
        % find daily return rate, then find mean daily return
        ret = [ret (stock_close(j+1,i) - stock_close(j,i))/stock_close(j,2)];
    end
    % yearly return rate is daily return rate compounded over a year
    yearly_ret = [yearly_ret (1+mean(ret))^(360)];
    retrec = [retrec; ret];
    daily_ret = [daily_ret mean(ret)];
end
%% Find Covariance Matrix
CovMat = zeros(size((stock_close),2),size((stock_close),2));
for i = 1:size((stock_close),2)
    for j = 1:size((stock_close),2)
        for k = 1:size((stock_close),1)-1
            CovMat(i,j) = CovMat(i,j) + (retrec(i,k)-daily_ret(i))*(retrec(j,k)-daily_ret(j));
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
n_hist = 100; % number of histogram bins
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
